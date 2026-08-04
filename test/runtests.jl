#!/usr/bin/env julia

# Regression tests for the physics invariants this project relies on.
#
#   julia --project=. test/runtests.jl
#   julia -t auto --project=. test/runtests.jl     (also exercises the threaded KPM)
#
# These are not unit tests of every function. They pin down the specific results
# that were validated by hand during development, so that a refactor cannot break
# them silently. Each one has been checked against an independent source: an
# analytical formula, a published table, or a mathematical identity.
#
# Deliberately excluded: anything that needs the experimental data files or takes
# more than a few minutes, so this stays runnable as a pre-commit check.

using Test
using Printf
using LinearAlgebra
using Statistics
using Random
using Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"));       using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

controls() = SV.sv_load_controls(REPO_ROOT)
params() = SV.sv_load_params(REPO_ROOT, controls()).params

@testset "YbZn2GaO5" begin

# -----------------------------------------------------------------------------
@testset "exchange shell geometry" begin
    # The crystal is P1, so Sunny will not generate the triangular shells from one
    # representative bond; the offsets are hard-coded. These sums must reproduce
    # the analytical J(q) form factors, or every dispersion is wrong.
    for row in SV.sv_exchange_geometry_sanity_rows()
        @test row.delta1_analytical ≈ row.delta1_bonds atol = 1e-10
        @test row.delta2_analytical ≈ row.delta2_bonds atol = 1e-10
    end
    # Known values at the high-symmetry points.
    K = [1 / 3, 1 / 3, 0.0]
    M = [0.5, 0.0, 0.0]
    @test SV.sv_analytical_delta1(K) ≈ 9.0 atol = 1e-10
    @test SV.sv_analytical_delta2(K) ≈ 0.0 atol = 1e-10
    @test SV.sv_analytical_delta1(M) ≈ 8.0 atol = 1e-10
    @test SV.sv_analytical_delta2(M) ≈ 8.0 atol = 1e-10
end

# -----------------------------------------------------------------------------
@testset "supercell construction" begin
    c = controls()
    p = params()
    # seed_dims must tile the cell, or fall back to a direct build.
    @test SV.sv_seed_divides_cell((36, 36, 1), (3, 3, 1))
    @test !SV.sv_seed_divides_cell((16, 16, 1), (3, 3, 1))
    @test SV.sv_repeat_factor((36, 36, 1), (3, 3, 1)) == (12, 12, 1)

    # A cell incommensurate with the seed must still build, via the direct path.
    b = SV.sv_build_supercell_system(p, c; cell_size=(8, 8, 1), seed_dims=(3, 3, 1))
    @test b.built_directly
    @test b.sys.dims == (8, 8, 1)

    b2 = SV.sv_build_supercell_system(p, c; cell_size=(6, 6, 1), seed_dims=(3, 3, 1))
    @test !b2.built_directly
    @test b2.sys.dims == (6, 6, 1)

    # realization = 0 must reproduce the original single-realization seed exactly,
    # or previously published neutron/KPM results silently change.
    g(b) = [b.sys.gs[s][3, 3] for s in eachsite(b.sys)]
    a0 = SV.sv_build_supercell_system(p, c; cell_size=(6, 6, 1), realization=0)
    a0b = SV.sv_build_supercell_system(p, c; cell_size=(6, 6, 1), realization=0)
    a1 = SV.sv_build_supercell_system(p, c; cell_size=(6, 6, 1), realization=1)
    @test g(a0) == g(a0b)          # deterministic
    @test g(a0) != g(a1)           # and realizations really differ
end

# -----------------------------------------------------------------------------
@testset "saturated moment" begin
    # Far above saturation every spin aligns, so the moment per site must equal
    # mean(g_i) * S for that realization. This ties together the g-tensor
    # convention, moment_sign, and sv_m_parallel_uB_per_site.
    c = controls()
    p = params()
    uhat = SV.sv_field_direction(c)
    spin_S = Float64(c["common"]["spin_S"])
    b = SV.sv_build_supercell_system(p, c; cell_size=(6, 6, 1), field_T=60.0)
    sys = b.sys
    SV.sv_set_field_T!(sys, uhat, b.units, 60.0)
    SV.sv_polarize_along_field!(sys, uhat; field_T=60.0)
    minimize_energy!(sys; maxiters=20_000)
    M = -SV.sv_m_parallel_uB_per_site(sys, uhat)     # moment_sign = -1
    Msat = SV.sv_m_sat_uB_per_site(sys, uhat, spin_S)
    @test M ≈ Msat rtol = 1e-4
    @test M > 0                                       # sign convention
end

# -----------------------------------------------------------------------------
@testset "two-component nonnegative scale fit" begin
    # Exactness on a constructed problem, and the nonnegativity guarantee.
    B = collect(range(0.5, 7.0; length=25))
    x = @. 1 - exp(-B / 2)
    a_true, b_true = 0.7, 0.05
    a, b = SV.sv_best_two_component_scale(a_true .* x .+ b_true .* B, x, B)
    @test a ≈ a_true rtol = 1e-8
    @test b ≈ b_true rtol = 1e-8
    # A target with no linear content must not produce a negative slope.
    a2, b2 = SV.sv_best_two_component_scale(x, x, B)
    @test a2 >= 0 && b2 >= 0
    # NaNs in the target must be ignored, not poison the fit.
    y = a_true .* x .+ b_true .* B
    y[[3, 11]] .= NaN
    a3, b3 = SV.sv_best_two_component_scale(y, x, B)
    @test a3 ≈ a_true rtol = 1e-6
    @test b3 ≈ b_true rtol = 1e-6
end

# -----------------------------------------------------------------------------
@testset "energy resolution kernel headroom" begin
    # The bug this guards: if the KPM kernel is wider than the target resolution,
    # the quadrature subtraction has nothing to remove and the model stays
    # over-broadened. That used to be clamped silently.
    c = controls()
    er = SV.sv_energy_resolution_controls(c; section="kpm")
    if er.enabled && er.subtract_kpm_kernel
        # Check the FULL configured energy grid, not just the fit window: the target
        # FWHM is smallest at high energy transfer, so the window alone would miss a
        # kernel that is too wide at the top of the range.
        kc = c["kpm"]
        egrid = range(Float64(kc["energy_min_meV"]), Float64(kc["energy_max_meV"]);
                      length=Int(kc["n_energy"]))
        h = SV.sv_resolution_kernel_headroom(er; energies=egrid)
        @test h.ok
        @test er.kernel_fwhm <= h.max_valid_kernel_fwhm + 1e-12
        # And the window on its own must also be satisfied.
        @test SV.sv_resolution_kernel_headroom(er; energies=range(0.5, 3.0; length=26)).ok
        # And the detector must actually fire on a deliberately too-wide kernel.
        bad = merge(er, (; kernel_fwhm=10.0))
        @test !SV.sv_resolution_kernel_headroom(bad; energies=range(0.5, 3.0; length=26)).ok
    end
end

# -----------------------------------------------------------------------------
@testset "KPM regularization for disordered systems" begin
    # Disorder puts magnon modes near zero energy, so Sunny's default
    # regularization of 1e-8 loses positive-definiteness and KPM aborts with
    # "Not an energy-minimum". This is NOT a relaxation failure: the state below
    # converges and its energy is correct. If someone drops [kpm].regularization
    # back to the Sunny default, this catches it.
    c = controls()
    p = merge(params(), (; J1_meV=0.25, J2_meV=0.01, sigma_J=0.5, gzz=3.8, sigma_gzz=0.8))
    uhat = SV.sv_field_direction(c)
    b = SV.sv_build_supercell_system(p, c; cell_size=(12, 12, 1), seed_dims=(3, 3, 1),
                                     field_T=9.0, realization=0)
    sys = b.sys
    SV.sv_set_field_T!(sys, uhat, b.units, 9.0)
    SV.sv_polarize_along_field!(sys, uhat; field_T=9.0)
    res = minimize_energy!(sys; maxiters=20_000)
    @test occursin("Converged", string(res))          # relaxation is not the problem

    E = collect(range(0.0, 4.0; length=21))
    qs = [[0.44, -0.12, 0.0], [0.5, 0.0, 0.0]]
    runs(reg) = try
        swt = SpinWaveTheoryKPM(sys; measure=SV.sv_sunny_measure(sys, c),
                                tol=0.05, regularization=reg)
        intensities(swt, qs; energies=E, kernel=gaussian(fwhm=0.05))
        true
    catch
        false
    end
    @test !runs(1e-8)                                  # the trap
    @test runs(1e-6)                                   # the fix
    @test Float64(get(c["kpm"], "regularization", 1e-8)) >= 1e-6
end

# -----------------------------------------------------------------------------
@testset "resolution deposition conserves weight at grid edges" begin
    # The bug: the old implementation renormalized the kernel over the target grid
    # unconditionally, so a model energy OUTSIDE the target range still deposited its
    # full intensity, piling onto the edge bins. That produced a spurious spike at the
    # top of the (0,1,0) cut, which ends at 3.275 meV while the model grid runs to
    # 4 meV. Interior bins must be unaffected; edge bins must not over-collect.
    c = controls()
    er = SV.sv_energy_resolution_controls(c; section="kpm")
    if er.enabled
        Em = collect(range(0.0, 4.0; length=161))
        Et = collect(range(0.025, 3.275; length=66))
        I = ones(length(Em), 1)
        out = SV.sv_post_deposit_energy_resolution(Em, I, Et, er)
        # A flat unit model over a uniform target grid deposits ~1 per model point per
        # target bin width; with 2 model points per target bin the interior is ~2.
        interior = 10:56
        @test all(1.5 .< out[interior, 1] .< 2.5)
        # The last bin must not be an order of magnitude above the interior.
        @test out[end, 1] < 2 * maximum(out[interior, 1])
        # Weight from far outside the range is dropped, not smeared in.
        lost = SV.sv_deposit_lost_fraction(Em, Et, er)
        @test lost[findfirst(>=(1.0), Em)] < 1e-6      # deep interior: nothing lost
        @test lost[findfirst(>=(3.8), Em)] > 0.99      # outside: all lost
        # Bin edges bracket the centres.
        edges = SV.sv_bin_edges_from_centers(Et)
        @test length(edges) == length(Et) + 1
        @test edges[1] < Et[1] && edges[end] > Et[end]
        @test issorted(edges)
    end
end

# -----------------------------------------------------------------------------
@testset "crystal field reproduces the published table" begin
    # Zhao et al., Phys. Rev. B 113, 014437 (2026), Tables II-IV. If Sunny's
    # Stevens normalization ever changed, this would catch it.
    B = Dict((2, 0) => -0.78, (4, 0) => 1.40e-2, (4, 3) => -8.2e-1,
             (6, 0) => 6.6e-4, (6, 3) => -3.0e-2, (6, 6) => 1.62e-2)
    O = stevens_matrices(7 / 2)
    H = Hermitian(sum(b .* O[n, m] for ((n, m), b) in B))
    E = sort(real.(eigvals(H)))
    E .-= E[1]
    levels = [mean(E[1:2]), mean(E[3:4]), mean(E[5:6]), mean(E[7:8])]
    @test levels[2] ≈ 38.3 atol = 1.5     # measured 38.3(1) meV
    @test levels[3] ≈ 60.6 atol = 1.5     # measured 60.6(1)
    @test levels[4] ≈ 95.4 atol = 1.5     # measured 95.4(2)
    # Ground-doublet g_par from the same eigenvectors, published as 3.44.
    F = eigen(H)
    V = F.vectors[:, sortperm(real.(F.values))]
    Jx, _, Jz = spin_matrices(7 / 2)
    blk = Hermitian(V[:, 1:2]' * Jz * V[:, 1:2])
    fb = eigen(blk)
    Vg = V[:, 1:2] * fb.vectors
    up = argmax(real.(fb.values))
    g_par = 2 * (8 / 7) * abs(real(Vg[:, up]' * Jz * Vg[:, up]))
    @test g_par ≈ 3.44 atol = 0.05
end

# -----------------------------------------------------------------------------
@testset "KPM threading is exact" begin
    # Threading over q must be bit-identical to serial, since q points are
    # independent and SpinWaveTheory clones the system.
    if Threads.nthreads() > 1
        c = controls()
        p = params()
        uhat = SV.sv_field_direction(c)
        b = SV.sv_build_supercell_system(p, c; cell_size=(6, 6, 1), field_T=14.0)
        sys = b.sys
        SV.sv_set_field_T!(sys, uhat, b.units, 14.0)
        SV.sv_polarize_along_field!(sys, uhat; field_T=14.0)
        minimize_energy!(sys; maxiters=20_000)

        energies = collect(range(0.0, 4.0; length=41))
        kern = gaussian(fwhm=0.08)
        meas() = SV.sv_sunny_measure(sys, c)
        rng = MersenneTwister(11)
        qs = [[0.5 + 0.05randn(rng), 0.05randn(rng), 0.0] for _ in 1:8]

        ser = intensities(SpinWaveTheoryKPM(sys; measure=meas(), tol=0.05,
                                            regularization=1e-5), qs;
                          energies, kernel=kern).data
        halves = [qs[1:2:end], qs[2:2:end]]
        outs = Vector{Any}(undef, 2)
        Threads.@threads for i in 1:2
            outs[i] = intensities(SpinWaveTheoryKPM(sys; measure=meas(), tol=0.05,
                                                   regularization=1e-5),
                                  halves[i]; energies, kernel=kern).data
        end
        @test Array(outs[1]) == Array(ser)[:, 1:2:end]
        @test Array(outs[2]) == Array(ser)[:, 2:2:end]
    else
        @info "KPM threading test skipped (needs julia -t auto)"
    end
end

# -----------------------------------------------------------------------------
@testset "configs parse and key invariants hold" begin
    cfgdir = joinpath(REPO_ROOT, "configs")
    for f in filter(x -> endswith(x, ".toml"), readdir(cfgdir; join=true))
        @test (load_toml_config(f); true)
    end
    c = controls()
    @test Float64(c["common"]["spin_S"]) == 0.5
    @test length(c["common"]["field_direction"]) == 3
    # Cell sizes used by the M(H) workflow must stay commensurate with the
    # three-sublattice 120-degree order, or the boundary frustrates it.
    for f in ("mvh_landscape_controls.toml", "mvh_fit_controls.toml")
        cfg = load_toml_config(joinpath(cfgdir, f))
        L = SV.sv_tuple3(cfg["run"]["cell_size"])[1]
        @test L % 3 == 0
    end
end

@testset "momentum resolution quadrature is exact" begin
    # The legacy rule placed n equally spaced nodes over +/- grid_nsigma*sigma and
    # renormalised Gaussian weights. Renormalising after truncation fixes the ZEROTH
    # moment but not the SECOND, so the realised resolution width was wrong -- 5.9% too
    # narrow at the production setting (n = 3, grid_nsigma = 1.5). Because momentum
    # resolution converts into energy width wherever the dispersion is steep, that biases
    # sigma_J high at K and M, which is precisely where sigma_J is determined.
    sigeff(xs, ws) = sqrt(sum(ws .* xs .^ 2) / sum(ws))

    # Gauss-Hermite must be exact for every node count, and must ignore grid_nsigma.
    for n in (2, 3, 5, 7, 9), nsig in (1.5, 3.0)
        xs, ws = SV.sv_gaussian_grid_axis(n, 2.5, nsig; quadrature=:gauss_hermite)
        @test isapprox(sum(ws), 1.0; atol=1e-12)
        @test isapprox(sigeff(xs, ws), 2.5; rtol=1e-10)
    end

    # n = 3 must be exactly {0, +/-sqrt(3)sigma} with weights {2/3, 1/6, 1/6}.
    xs, ws = SV.sv_gaussian_grid_axis(3, 1.0, 1.5)
    @test isapprox(sort(xs), [-sqrt(3), 0.0, sqrt(3)]; rtol=1e-12)
    @test isapprox(sort(ws), [1/6, 1/6, 2/3]; rtol=1e-12)

    # Pin the legacy rule's error so nobody "fixes" the default back to it. Both obvious
    # repairs make it worse: widening the window at n = 3 is catastrophic because three
    # nodes at +/-3 sigma put ~98% of the weight at the centre, and adding nodes at
    # grid_nsigma = 1.5 crowds them into a window that is already too narrow.
    for (n, nsig, expected) in ((3, 1.5, 0.941164), (3, 3.0, 0.442285),
                                (5, 1.5, 0.855154), (9, 1.5, 0.802254))
        xs, ws = SV.sv_gaussian_grid_axis(n, 1.0, nsig; quadrature=:truncated_gaussian_grid)
        @test isapprox(sigeff(xs, ws), expected; rtol=1e-4)
    end

    # The default must be the exact rule, on both the 1D and 2D control paths.
    c1 = SV.sv_kpm_1d_q_averaging_controls(Dict("kpm" => Dict("q_averaging" => Dict())))
    @test c1.quadrature === :gauss_hermite
    c2 = SV.sv_kpm_2d_q_averaging_controls(Dict("q_averaging" => Dict()))
    @test c2.quadrature === :gauss_hermite
    @test SV.sv_kpm_1d_q_averaging_controls(Dict("kpm" => Dict("q_averaging" =>
        Dict("resolution_quadrature" => "truncated_gaussian_grid")))).quadrature ===
        :truncated_gaussian_grid
end

@testset "cached KPM operators are bit-identical to fresh builds" begin
    # SpinWaveTheoryKPM construction is ~0.36 s at 36x36x1 and was happening once per
    # chunk on EVERY intensities call, so a context serving several cuts rebuilt identical
    # operators for each. Reuse is only sound if the object carries no state across calls;
    # the plausible failure is buffers sized by the first call's q count, which would
    # silently corrupt later calls. Bit-identity is the acceptance criterion -- anything
    # less and the speedup is not worth having.
    controls = SV.sv_load_controls(REPO_ROOT)
    controls["kpm"]["system_size"] = [12, 12, 1]
    controls["kpm"]["dims"] = [3, 3, 1]
    controls["kpm"]["repeat_factor"] = [4, 4, 1]
    controls["kpm"]["maxiters"] = 500
    controls["kpm"]["regularization"] = 1e-5
    (; params) = SV.sv_load_params(REPO_ROOT, controls)
    params = merge(params, (; J1_meV=0.25, J2_meV=0.01, sigma_J=0.5,
                              gzz=3.35, sigma_gzz=0.8))
    ctx = SV.sv_kpm_context(params, controls; component=:dispersive, field_T=9.0,
                            realization=0, maxiters=500, relax_attempts=1)
    @test !isnothing(get(ctx, :swt_cache, nothing))   # on by default

    energies = collect(range(0.0, 4.0; length=21))
    kern = gaussian(fwhm=0.05)
    mkqs(n) = [[0.33 + 0.01*i, 0.33 + 0.005*i, 0.0] for i in 1:n]
    function fill_it(qs, nchunks, cache)
        I0 = zeros(Float64, length(energies), length(qs))
        SV._sv_kpm_fill_intensity!(I0, ctx.sys, controls, qs, energies, kern, 0.05, 1e-5,
                                   nchunks, length(energies), length(qs); cache)
        return I0
    end

    # Repeated calls through one cached pool must match fresh builds exactly.
    cache = Dict{Any,Any}()
    qs = mkqs(12)
    ref = fill_it(qs, 4, nothing)
    for _ in 1:3
        @test fill_it(qs, 4, cache) == ref
    end

    # The dangerous case: DIFFERENT q counts reusing the same pooled operators.
    cache2 = Dict{Any,Any}()
    for n in (12, 6, 18, 12)
        qsn = mkqs(n)
        @test fill_it(qsn, 4, cache2) == fill_it(qsn, 4, nothing)
    end

    # Keyed on nchunks, so changing the chunk count must build a new pool rather than
    # reuse a wrongly sized one. Three distinct chunk counts => exactly three keys.
    cache3 = Dict{Any,Any}()
    for nch in (1, 4, 8, 4, 1)
        @test fill_it(qs, nch, cache3) == fill_it(qs, nch, nothing)
    end
    @test length(cache3) == 3

    # The escape hatch must actually disable it.
    controls["kpm"]["cache_kpm_operators"] = false
    ctx2 = SV.sv_kpm_context(params, controls; component=:dispersive, field_T=9.0,
                             realization=0, maxiters=500, relax_attempts=1)
    @test isnothing(get(ctx2, :swt_cache, nothing))
end

end
