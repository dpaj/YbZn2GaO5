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

        ser = intensities(SpinWaveTheoryKPM(sys; measure=meas(), tol=0.05), qs;
                          energies, kernel=kern).data
        halves = [qs[1:2:end], qs[2:2:end]]
        outs = Vector{Any}(undef, 2)
        Threads.@threads for i in 1:2
            outs[i] = intensities(SpinWaveTheoryKPM(sys; measure=meas(), tol=0.05),
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

end
