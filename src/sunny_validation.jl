module SunnyValidation

using LinearAlgebra
using Random
using Printf
using StaticArrays
using DelimitedFiles
using Statistics

using Sunny
using CairoMakie

using Main.YZGOCofit

const SV_MU_B_MEV_PER_T = 0.05788381806
const SV_KB_MEV_PER_K = 0.08617333262

# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------

sv_symbol(x::AbstractString) = Symbol(x)
sv_symbol(x::Symbol) = x

function sv_repo_path(repo_root::AbstractString, p::AbstractString)
    isabspath(p) && return normpath(p)
    return normpath(joinpath(repo_root, splitpath(p)...))
end

function sv_load_controls(repo_root::AbstractString)
    path = joinpath(repo_root, "configs", "sunny_validation_controls.toml")
    return load_toml_config(path)
end

function sv_load_params(repo_root::AbstractString, controls::Dict)
    pth = sv_repo_path(repo_root, controls["paths"]["best_fit_parameters_toml"])
    params = load_canonical_model_parameters(pth)
    return (; params, path=pth)
end

function sv_field_grid(controls::Dict)
    c = controls["magnetization"]
    return collect(range(Float64(c["Bmin_T"]), Float64(c["Bmax_T"]); length=Int(c["nB"])))
end

function sv_second_kernel_weight(params, controls::Dict)
    if get(controls["weights"], "use_second_kernel_relative_intensity_from_best_fit", true)
        return params.second_kernel_relative_intensity
    else
        return Float64(controls["weights"]["manual_second_kernel_relative_intensity"])
    end
end

function sv_field_direction(controls::Dict)
    u = Float64.(controls["common"]["field_direction"])
    n = norm(u)
    n > 0 || error("field_direction must be nonzero")
    return SVector{3,Float64}(u ./ n)
end

function sv_dims_tuple(v)
    length(v) == 3 || error("Expected a length-3 dims/repeat_factor vector")
    return (Int(v[1]), Int(v[2]), Int(v[3]))
end

function sv_units()
    return Units(:meV, :angstrom)
end

function sv_try_pkgversion(mod)
    try
        return string(pkgversion(mod))
    catch
        return "unknown"
    end
end

function sv_read_magnetization_csv(path::AbstractString)
    B = Float64[]
    M = Float64[]
    isfile(path) || error("Could not find magnetization CSV: $path")
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            parts = split(s, ',')
            length(parts) < 2 && continue
            b = tryparse(Float64, strip(parts[1]))
            m = tryparse(Float64, strip(parts[2]))
            if b !== nothing && m !== nothing
                push!(B, b)
                push!(M, m)
            end
        end
    end
    isempty(B) && error("No numeric magnetization data in $path")
    idx = sortperm(B)
    return (; B_T=B[idx], M_muB_per_Yb=M[idx], path=path)
end

function sv_interp1(x::AbstractVector, y::AbstractVector, xq::AbstractVector)
    out = fill(NaN, length(xq))
    length(x) == length(y) || error("x and y length mismatch")
    for (i, q) in enumerate(xq)
        if q < first(x) || q > last(x)
            continue
        end
        j = searchsortedlast(x, q)
        if j <= 0
            out[i] = first(y)
        elseif j >= length(x)
            out[i] = last(y)
        else
            t = (q - x[j]) / (x[j+1] - x[j])
            out[i] = (1-t)*y[j] + t*y[j+1]
        end
    end
    return out
end

function sv_best_vertical_scale(model::AbstractVector, data::AbstractVector)
    mask = isfinite.(model) .& isfinite.(data)
    count(mask) == 0 && return NaN
    den = sum(abs2, model[mask])
    den <= eps(Float64) && return NaN
    return sum(model[mask] .* data[mask]) / den
end


# -----------------------------------------------------------------------------
# Magnetization scale and two-component mixture convention
# -----------------------------------------------------------------------------

function sv_magnetization_global_scale(params, controls::Dict)
    mag = controls["magnetization"]
    if get(mag, "use_magnetization_global_scale_from_best_fit", true)
        if hasproperty(params, :magnetization_global_scale)
            return Float64(params.magnetization_global_scale)
        else
            error(
                "best_fit_parameters.toml does not define magnetization_global_scale. " *
                "Add [magnetization_extrinsic] magnetization_global_scale = <value>."
            )
        end
    else
        return Float64(get(mag, "manual_magnetization_global_scale", 1.0))
    end
end

function sv_magnetization_mixture(Mdisp::AbstractVector, Mflat::AbstractVector, r2::Real, controls::Dict)
    normalize_weight = get(controls["magnetization"], "normalize_second_kernel_weight", true)
    if normalize_weight
        denom = 1.0 + Float64(r2)
        return (Mdisp ./ denom, (Float64(r2) .* Mflat) ./ denom)
    else
        return (copy(Mdisp), Float64(r2) .* Mflat)
    end
end

function sv_build_magnetization_comparison(Bs, Mdisp_raw, Mflat_raw, params, controls::Dict; moment_sign::Real=1.0)
    sign_factor = Float64(moment_sign)

    # Convert calculator convention to the experimental positive-M convention.
    Mdisp_site = sign_factor .* Mdisp_raw
    Mflat_site = sign_factor .* Mflat_raw

    # Mixture convention mirrors the analytical co-fit comparison layer.
    r2 = sv_second_kernel_weight(params, controls)
    Mdisp_contrib_unscaled, Mflat_contrib_unscaled =
        sv_magnetization_mixture(Mdisp_site, Mflat_site, r2, controls)

    Mmag_unscaled = Mdisp_contrib_unscaled .+ Mflat_contrib_unscaled

    include_chi_vv = get(controls["magnetization"], "include_chi_vv", true)
    Mvv_unscaled = include_chi_vv ?
        params.chi_vv_muB_per_T .* Bs :
        zeros(Float64, length(Bs))

    Mcombo_unscaled = Mmag_unscaled .+ Mvv_unscaled

    # Analytical-model convention:
    # magnetization_global_scale multiplies the whole unscaled magnetization
    # comparison curve, including the Van Vleck term.
    magnetization_scale = sv_magnetization_global_scale(params, controls)
    Mtotal_scaled = magnetization_scale .* Mcombo_unscaled

    return (;
        r2 = r2,
        magnetization_scale = magnetization_scale,
        normalize_second_kernel_weight = get(controls["magnetization"], "normalize_second_kernel_weight", true),
        include_chi_vv = include_chi_vv,
        moment_sign = sign_factor,

        M_disp_raw_uB_per_site = Mdisp_raw,
        M_flat_raw_uB_per_site = Mflat_raw,
        M_disp_signed_uB_per_site = Mdisp_site,
        M_flat_signed_uB_per_site = Mflat_site,

        M_disp_contrib_unscaled = Mdisp_contrib_unscaled,
        M_flat_contrib_unscaled = Mflat_contrib_unscaled,
        M_magnetic_unscaled = Mmag_unscaled,
        M_vv_unscaled = Mvv_unscaled,
        M_combo_unscaled = Mcombo_unscaled,

        M_disp_scaled = magnetization_scale .* Mdisp_contrib_unscaled,
        M_flat_scaled = magnetization_scale .* Mflat_contrib_unscaled,
        M_magnetic_scaled = magnetization_scale .* Mmag_unscaled,
        M_vv_scaled = magnetization_scale .* Mvv_unscaled,
        M_total_scaled = Mtotal_scaled,
    )
end

function sv_magnetization_scale_diagnostics(Bs, combo_unscaled, total_scaled, data)
    data_on_grid = sv_interp1(data.B_T, data.M_muB_per_Yb, Bs)
    return (;
        diagnostic_scale_unscaled_combo = sv_best_vertical_scale(combo_unscaled, data_on_grid),
        diagnostic_scale_primary_total = sv_best_vertical_scale(total_scaled, data_on_grid),
        M_exp_interp_uB_per_Yb = data_on_grid,
    )
end

function sv_write_magnetization_csv(path::AbstractString, Bs, comp, diag)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join([
            "B_T",
            "M_total_uB_per_Yb",
            "M_magnetic_scaled",
            "M_disp_scaled",
            "M_flat_scaled",
            "M_vv_scaled",
            "M_combo_unscaled",
            "M_magnetic_unscaled",
            "M_disp_contrib_unscaled",
            "M_flat_contrib_unscaled",
            "M_vv_unscaled",
            "M_disp_raw_uB_per_site",
            "M_flat_raw_uB_per_site",
            "M_exp_interp_uB_per_Yb",
            "second_kernel_weight",
            "magnetization_global_scale",
            "diagnostic_scale_unscaled_combo",
            "diagnostic_scale_primary_total",
            "normalize_second_kernel_weight",
            "include_chi_vv",
            "moment_sign",
        ], ","))

        for i in eachindex(Bs)
            vals = [
                Bs[i],
                comp.M_total_scaled[i],
                comp.M_magnetic_scaled[i],
                comp.M_disp_scaled[i],
                comp.M_flat_scaled[i],
                comp.M_vv_scaled[i],
                comp.M_combo_unscaled[i],
                comp.M_magnetic_unscaled[i],
                comp.M_disp_contrib_unscaled[i],
                comp.M_flat_contrib_unscaled[i],
                comp.M_vv_unscaled[i],
                comp.M_disp_raw_uB_per_site[i],
                comp.M_flat_raw_uB_per_site[i],
                diag.M_exp_interp_uB_per_Yb[i],
                comp.r2,
                comp.magnetization_scale,
                diag.diagnostic_scale_unscaled_combo,
                diag.diagnostic_scale_primary_total,
                comp.normalize_second_kernel_weight ? 1.0 : 0.0,
                comp.include_chi_vv ? 1.0 : 0.0,
                comp.moment_sign,
            ]
            println(io, join([@sprintf("%.10g", Float64(v)) for v in vals], ","))
        end
    end
    return path
end

function sv_write_xy_csv(path::AbstractString, header::AbstractString, cols...)
    mkpath(dirname(path))
    n = length(cols[1])
    open(path, "w") do io
        println(io, header)
        for i in 1:n
            for (j, c) in enumerate(cols)
                j > 1 && print(io, ',')
                print(io, @sprintf("%.10g", Float64(c[i])))
            end
            println(io)
        end
    end
    return path
end

# -----------------------------------------------------------------------------
# Current analytical ingredients reused for validation targets
# -----------------------------------------------------------------------------

function sv_exchange_Dmax_meV(J1_meV::Real, J2_meV::Real; mode::Symbol=:high_symmetry)
    J1 = Float64(J1_meV)
    J2 = Float64(J2_meV)
    if mode == :high_symmetry
        return max(0.0, 9.0 * J1, 8.0 * (J1 + J2))
    elseif mode == :K_only
        return max(0.0, 9.0 * J1)
    elseif mode == :M_only
        return max(0.0, 8.0 * (J1 + J2))
    else
        error("Unknown Dmax mode $mode")
    end
end

function sv_normal_quadrature(n::Int=301; zmax::Real=5.0)
    n >= 5 || error("normal quadrature needs at least 5 points")
    z = collect(range(-Float64(zmax), Float64(zmax); length=n))
    w = exp.(-0.5 .* z .^ 2)
    w ./= sum(w)
    return (; z, w)
end

function sv_make_draws(n::Int; seed::Int=20260611)
    rng = MersenneTwister(seed)
    return (; zg=randn(rng, n), z1=randn(rng, n), z2=randn(rng, n))
end

function sv_dispersive_magnetization_analytic(Bs_T, params, draws; S::Real=0.5)
    Sval = Float64(S)
    n = length(draws.zg)
    Bsat = Vector{Float64}(undef, n)
    msat = Vector{Float64}(undef, n)

    sJ1_abs = params.sigma_J * abs(params.J1_meV)
    sJ2_abs = params.sigma_J * abs(params.J2_meV)
    for i in 1:n
        gi = params.gzz + params.sigma_gzz * draws.zg[i]
        J1i = params.J1_meV + sJ1_abs * draws.z1[i]
        J2i = params.J2_meV + sJ2_abs * draws.z2[i]
        Dmaxi = sv_exchange_Dmax_meV(J1i, J2i)
        msat[i] = gi * Sval
        Bsat[i] = (isfinite(gi) && gi > 0 && Dmaxi > 0) ? Sval * Dmaxi / (gi * SV_MU_B_MEV_PER_T) : NaN
    end

    out = zeros(Float64, length(Bs_T))
    for (k, B) in enumerate(Bs_T)
        total = 0.0
        nvalid = 0
        for i in 1:n
            if isfinite(Bsat[i]) && isfinite(msat[i]) && Bsat[i] >= 0
                nvalid += 1
                if Bsat[i] <= 0
                    total += sign(B) * msat[i]
                else
                    total += sign(B) * msat[i] * min(abs(B)/Bsat[i], 1.0)
                end
            end
        end
        out[k] = nvalid > 0 ? total / nvalid : NaN
    end
    return out
end

function sv_flat_magnetization_independent_spin(Bs_T, params, quad; temperature_K::Real=0.42, S::Real=0.5)
    Bq = Float64.(Bs_T)
    T = Float64(temperature_K)
    Sval = Float64(S)
    gvals_all = params.gzz2 .+ params.sigma_gzz2 .* quad.z
    valid = isfinite.(gvals_all) .& (gvals_all .> 0.0)
    count(valid) > 0 || return fill(NaN, length(Bq))
    gvals = gvals_all[valid]
    w = quad.w[valid]
    w ./= sum(w)
    out = zeros(Float64, length(Bq))
    if T <= 0
        for (i, B) in enumerate(Bq)
            out[i] = sign(B) * sum(w .* (gvals .* Sval))
        end
    else
        denom = SV_KB_MEV_PER_K * T
        for (i, B) in enumerate(Bq)
            acc = 0.0
            for j in eachindex(gvals)
                arg = gvals[j] * SV_MU_B_MEV_PER_T * B * Sval / denom
                acc += w[j] * gvals[j] * Sval * tanh(arg)
            end
            out[i] = acc
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Mean-field bridge validation
# -----------------------------------------------------------------------------

function sv_run_meanfield_magnetization(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    print_canonical_model_parameters(params)
    @info "Sunny validation: mean-field bridge" sunny_version=sv_try_pkgversion(Sunny) params_path=path

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    out_fig_dir = sv_repo_path(repo_root, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir); mkpath(out_fig_dir)

    Bs = sv_field_grid(controls)
    mag_path = sv_repo_path(repo_root, controls["paths"]["magnetization_csv"])
    data = sv_read_magnetization_csv(mag_path)

    S = Float64(controls["common"]["spin_S"])
    seed = Int(controls["common"]["seed"])
    n_samples = Int(controls["meanfield"]["n_disorder_samples"])
    qn = Int(controls["meanfield"]["normal_quad_n"])
    qzmax = Float64(controls["meanfield"]["normal_quad_zmax"])
    T = Float64(controls["magnetization"]["temperature_K"])
    moment_sign = Float64(get(controls["magnetization"], "moment_sign", 1.0))

    draws = sv_make_draws(n_samples; seed)
    quad = sv_normal_quadrature(qn; zmax=qzmax)

    Mdisp_raw = sv_dispersive_magnetization_analytic(Bs, params, draws; S)
    Mflat_raw = sv_flat_magnetization_independent_spin(Bs, params, quad; temperature_K=T, S)

    comp = sv_build_magnetization_comparison(Bs, Mdisp_raw, Mflat_raw, params, controls; moment_sign)
    diag = sv_magnetization_scale_diagnostics(Bs, comp.M_combo_unscaled, comp.M_total_scaled, data)

    @info "Magnetization scale convention: mean-field"         magnetization_global_scale=comp.magnetization_scale         diagnostic_scale_unscaled_combo=diag.diagnostic_scale_unscaled_combo         diagnostic_scale_primary_total=diag.diagnostic_scale_primary_total         normalize_weight=comp.normalize_second_kernel_weight         moment_sign=comp.moment_sign

    csv_path = joinpath(out_table_dir, "sunny_meanfield_magnetization.csv")
    sv_write_magnetization_csv(csv_path, Bs, comp, diag)

    fig = Figure(size=(950, 600))
    ax = Axis(fig[1,1], xlabel="B (T)", ylabel="M (μB / Yb)", title="Sunny validation: mean-field bridge")
    scatter!(ax, data.B_T, data.M_muB_per_Yb, label="experiment")
    lines!(ax, Bs, comp.M_total_scaled, label="total × analytical scale")
    lines!(ax, Bs, comp.M_disp_scaled, label="dispersive contrib")
    lines!(ax, Bs, comp.M_flat_scaled, label="flat contrib")
    lines!(ax, Bs, comp.M_vv_scaled, label="Van Vleck contrib")

    if get(controls["magnetization"], "plot_diagnostic_scaled_curve", true)
        lines!(ax, Bs, diag.diagnostic_scale_unscaled_combo .* comp.M_combo_unscaled,
            label=@sprintf("diagnostic scaled total ×%.3g", diag.diagnostic_scale_unscaled_combo))
    end

    axislegend(ax, position=:rb)
    fig_path = joinpath(out_fig_dir, "sunny_meanfield_magnetization.png")
    save(fig_path, fig)

    @info "Saved mean-field Sunny validation" csv_path fig_path scale=comp.magnetization_scale normalize_weight=comp.normalize_second_kernel_weight moment_sign=comp.moment_sign
    return (; B_T=Bs, M_total=comp.M_total_scaled, M_disp=comp.M_disp_scaled,
        M_flat=comp.M_flat_scaled, M_vv=comp.M_vv_scaled, r2=comp.r2,
        scale=comp.magnetization_scale, diagnostics=diag, csv_path, fig_path)
end

# -----------------------------------------------------------------------------
# Preliminary effective triangular Sunny system builder
# -----------------------------------------------------------------------------

function sv_effective_triangle_crystal(controls::Dict)
    a = Float64(controls["common"]["lattice_a_angstrom"])
    c = Float64(controls["common"]["lattice_c_angstrom"])
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    return Crystal(latvecs, [[0.0, 0.0, 0.0]], 1; types=["Yb"])
end

function sv_set_gzz_all_sites!(sys, gzz::Real)
    for site in eachsite(sys)
        G = sys.gs[site]
        Gm = Matrix(G)
        Gm .= 0.0
        Gm[1,1] = 0.0
        Gm[2,2] = 0.0
        Gm[3,3] = Float64(gzz)
        sys.gs[site] = SMatrix{3,3,Float64,9}(Gm)
    end
    return sys
end

function sv_build_effective_sunny_system(params, controls::Dict; component::Symbol=:dispersive, dims=(4,4,1), field_T::Real=0.0)
    units = sv_units()
    cryst = sv_effective_triangle_crystal(controls)
    g0 = component == :flat ? params.gzz2 : params.gzz
    sys = System(cryst, [1 => Moment(s=Float64(controls["common"]["spin_S"]), g=Float64(g0))], :dipole; dims=dims, seed=Int(controls["common"]["seed"]))
    # Force Ising-like g convention for H || c validation.  This may need tuning
    # if Sunny's current Moment constructor already created an acceptable matrix.
    sv_set_gzz_all_sites!(sys, g0)

    if component == :dispersive
        # Effective triangular lattice.  These three directions generate the six
        # nearest neighbors with Hermitian counterparts.  J2 directions are a
        # first pass consistent with the analytical J1/J2 parameterization.
        for off in ([1,0,0], [0,1,0], [1,-1,0])
            set_exchange!(sys, params.J1_meV, Bond(1, 1, off))
        end
        for off in ([1,1,0], [2,-1,0], [1,-2,0])
            set_exchange!(sys, params.J2_meV, Bond(1, 1, off))
        end
    elseif component == :flat
        # Zero exchange: independent-spin/flat component.
    else
        error("Unknown component $component")
    end

    u = sv_field_direction(controls)
    set_field!(sys, collect(u .* (Float64(field_T) * units.T)))
    return (; sys, crystal=cryst, units)
end

function sv_apply_disorder!(sys, params, controls::Dict; component::Symbol=:dispersive, include_exchange::Bool=true, include_gzz::Bool=true)
    rng = MersenneTwister(Int(controls["common"]["seed"]) + 7919 + (component == :flat ? 17 : 0))
    if component == :dispersive && include_exchange
        for off in ([1,0,0], [0,1,0], [1,-1,0])
            for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
                set_exchange_at!(sys, params.J1_meV * (1 + params.sigma_J * randn(rng)), s1, s2; offset=o)
            end
        end
        for off in ([1,1,0], [2,-1,0], [1,-2,0])
            for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
                set_exchange_at!(sys, params.J2_meV * (1 + params.sigma_J * randn(rng)), s1, s2; offset=o)
            end
        end
    end
    if include_gzz
        sigma = component == :flat ? params.sigma_gzz2 : params.sigma_gzz
        gmean = component == :flat ? params.gzz2 : params.gzz
        for site in eachsite(sys)
            gi = max(0.0, gmean + sigma * randn(rng))
            G = sys.gs[site]
            Gm = Matrix(G)
            Gm[3,3] = gi
            sys.gs[site] = SMatrix{3,3,Float64,9}(Gm)
        end
    end
    return sys
end

function sv_m_parallel_uB_per_site(sys, uhat::SVector{3,Float64})
    acc = 0.0
    n = 0
    for site in eachsite(sys)
        m = sys.gs[site] * sys.dipoles[site]
        acc += dot(m, uhat)
        n += 1
    end
    return n > 0 ? acc / n : NaN
end

function sv_sweep_largecell_component(params, controls::Dict; component::Symbol, Bs_T, dims, repeat_factor, include_exchange_disorder::Bool, include_gzz_disorder::Bool, maxiters::Int)
    # Build at zero field, then continue field sweep.
    base = sv_build_effective_sunny_system(params, controls; component, dims, field_T=0.0)
    sys = base.sys
    if repeat_factor != (1,1,1)
        sys = to_inhomogeneous(repeat_periodically(sys, repeat_factor))
    else
        sys = to_inhomogeneous(sys)
    end
    sv_apply_disorder!(sys, params, controls; component, include_exchange=include_exchange_disorder, include_gzz=include_gzz_disorder)

    Random.seed!(Int(controls["common"]["seed"]) + 101)
    randomize_spins!(sys)
    minimize_energy!(sys; maxiters)

    units = sv_units()
    u = sv_field_direction(controls)
    M = zeros(Float64, length(Bs_T))
    for (i, B) in enumerate(Bs_T)
        set_field!(sys, collect(u .* (Float64(B) * units.T)))
        if get(controls["largecell"], "randomize_each_field", false)
            randomize_spins!(sys)
        end
        minimize_energy!(sys; maxiters)
        M[i] = sv_m_parallel_uB_per_site(sys, u)
        if i == 1 || i == length(Bs_T) || i % max(1, length(Bs_T) ÷ 10) == 0
            @printf("  %-10s B = %6.3f T -> M = % .8f μB/site\n", string(component), B, M[i])
        end
    end
    return M
end

function sv_run_largecell_magnetization(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    print_canonical_model_parameters(params)
    @info "Sunny validation: large-cell minimize_energy!" sunny_version=sv_try_pkgversion(Sunny) params_path=path

    Bs = sv_field_grid(controls)
    lc = controls["largecell"]
    dims = sv_dims_tuple(lc["dims"])
    repeat_factor = sv_dims_tuple(lc["repeat_factor"])
    include_exchange = get(lc, "include_exchange_disorder", true)
    include_gzz = get(lc, "include_gzz_disorder", true)
    maxiters = Int(lc["maxiters"])
    moment_sign = Float64(get(lc, "moment_sign", -1.0))

    Mdisp_raw = sv_sweep_largecell_component(params, controls; component=:dispersive, Bs_T=Bs, dims, repeat_factor, include_exchange_disorder=include_exchange, include_gzz_disorder=include_gzz, maxiters)
    Mflat_raw = sv_sweep_largecell_component(params, controls; component=:flat, Bs_T=Bs, dims, repeat_factor, include_exchange_disorder=false, include_gzz_disorder=include_gzz, maxiters)

    comp = sv_build_magnetization_comparison(Bs, Mdisp_raw, Mflat_raw, params, controls; moment_sign)

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    out_fig_dir = sv_repo_path(repo_root, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir); mkpath(out_fig_dir)

    mag_path = sv_repo_path(repo_root, controls["paths"]["magnetization_csv"])
    data = sv_read_magnetization_csv(mag_path)
    diag = sv_magnetization_scale_diagnostics(Bs, comp.M_combo_unscaled, comp.M_total_scaled, data)

    @info "Magnetization scale convention: large-cell"         magnetization_global_scale=comp.magnetization_scale         diagnostic_scale_unscaled_combo=diag.diagnostic_scale_unscaled_combo         diagnostic_scale_primary_total=diag.diagnostic_scale_primary_total         normalize_weight=comp.normalize_second_kernel_weight         moment_sign=comp.moment_sign

    csv_path = joinpath(out_table_dir, "sunny_largecell_magnetization.csv")
    sv_write_magnetization_csv(csv_path, Bs, comp, diag)

    fig = Figure(size=(950, 600))
    ax = Axis(fig[1,1], xlabel="B (T)", ylabel="M (μB / Yb)", title="Sunny validation: large-cell minimize_energy!")
    scatter!(ax, data.B_T, data.M_muB_per_Yb, label="experiment")
    lines!(ax, Bs, comp.M_total_scaled, label="total × analytical scale")
    lines!(ax, Bs, comp.M_disp_scaled, label="dispersive Sunny contrib")
    lines!(ax, Bs, comp.M_flat_scaled, label="flat Sunny contrib")
    lines!(ax, Bs, comp.M_vv_scaled, label="Van Vleck contrib")

    if get(controls["magnetization"], "plot_diagnostic_scaled_curve", true)
        lines!(ax, Bs, diag.diagnostic_scale_unscaled_combo .* comp.M_combo_unscaled,
            label=@sprintf("diagnostic scaled total ×%.3g", diag.diagnostic_scale_unscaled_combo))
    end

    axislegend(ax, position=:rb)
    fig_path = joinpath(out_fig_dir, "sunny_largecell_magnetization.png")
    save(fig_path, fig)

    @info "Saved large-cell Sunny validation" csv_path fig_path scale=comp.magnetization_scale normalize_weight=comp.normalize_second_kernel_weight moment_sign=comp.moment_sign
    return (; B_T=Bs, M_total=comp.M_total_scaled, M_disp=comp.M_disp_scaled,
        M_flat=comp.M_flat_scaled, M_vv=comp.M_vv_scaled, r2=comp.r2,
        scale=comp.magnetization_scale, diagnostics=diag, csv_path, fig_path)
end

# -----------------------------------------------------------------------------
# Preliminary KPM 1D spin-wave validation
# -----------------------------------------------------------------------------

function sv_qtag_to_q(qtag::AbstractString)
    if qtag == "0_1_0"
        return [0.0, 1.0, 0.0]
    elseif qtag == "0p33_0p33_0"
        return [1/3, 1/3, 0.0]
    elseif qtag == "0p5_0_0"
        return [0.5, 0.0, 0.0]
    else
        error("Unknown qtag $qtag")
    end
end

function sv_try_extract_sunny_intensity(res)
    # Sunny result internals can vary across versions.  Keep this deliberately
    # introspective for the first scaffold.  If this fails, inspect
    # propertynames(res) and patch this function.
    for nm in (:data, :intensity, :intensities, :I)
        if hasproperty(res, nm)
            return getproperty(res, nm)
        end
    end
    error("Could not extract intensity array from Sunny result. propertynames(res) = $(propertynames(res))")
end

function sv_kpm_component_spectrum(params, controls::Dict; component::Symbol, field_T::Real, qtag::AbstractString)
    kc = controls["kpm"]
    dims = sv_dims_tuple(kc["dims"])
    repeat_factor = sv_dims_tuple(kc["repeat_factor"])
    include_exchange = get(kc, "include_exchange_disorder", true)
    include_gzz = get(kc, "include_gzz_disorder", true)
    maxiters = Int(kc["maxiters"])

    base = sv_build_effective_sunny_system(params, controls; component, dims, field_T)
    sys = base.sys
    cryst = base.crystal
    if repeat_factor != (1,1,1)
        sys = to_inhomogeneous(repeat_periodically(sys, repeat_factor))
    else
        sys = to_inhomogeneous(sys)
    end
    sv_apply_disorder!(sys, params, controls; component, include_exchange=include_exchange, include_gzz=include_gzz)

    randomize_spins!(sys)
    minimize_energy!(sys; maxiters)

    q = sv_qtag_to_q(qtag)
    qs = [q]
    energies = collect(range(Float64(kc["energy_min_meV"]), Float64(kc["energy_max_meV"]); length=Int(kc["n_energy"])))
    kernel = gaussian(fwhm=Float64(kc["kernel_fwhm_meV"]))
    swt = SpinWaveTheoryKPM(sys; measure=ssf_perp(sys), tol=Float64(kc["tol"]))
    # Sunny supports either QPath objects or explicit arrays of q-points.
    # Use a single explicit q here because the experimental validation target is
    # a fixed 1D energy cut at each qtag, not a high-symmetry path plot.
    res = intensities(swt, qs; energies, kernel)
    raw = sv_try_extract_sunny_intensity(res)
    # For a fixed-q two-point path, average all q columns if possible.
    I = if raw isa AbstractVector
        Float64.(raw)
    elseif raw isa AbstractMatrix
        vec(mean(Float64.(raw); dims=1))
    else
        error("Unsupported Sunny intensity container type: $(typeof(raw))")
    end
    # If orientation is transposed, trim/interpolate by length.
    if length(I) != length(energies)
        I = vec(Float64.(raw))[1:min(end, length(energies))]
        if length(I) < length(energies)
            append!(I, fill(NaN, length(energies)-length(I)))
        end
    end
    return (; energy_meV=energies, intensity=I, result=res)
end

function sv_run_kpm_1d(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    print_canonical_model_parameters(params)
    @info "Sunny validation: KPM 1D spin-wave" sunny_version=sv_try_pkgversion(Sunny) params_path=path

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    out_fig_dir = sv_repo_path(repo_root, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir); mkpath(out_fig_dir)

    fields = Float64.(controls["common"]["fields_T"])
    qtags = String.(controls["kpm"]["qtags"])
    r2 = sv_second_kernel_weight(params, controls)

    all_rows = NamedTuple[]
    for B in fields
        for qtag in qtags
            @info "Computing Sunny KPM cut" B_T=B qtag=qtag
            disp = sv_kpm_component_spectrum(params, controls; component=:dispersive, field_T=B, qtag)
            flat = sv_kpm_component_spectrum(params, controls; component=:flat, field_T=B, qtag)
            total = disp.intensity .+ r2 .* flat.intensity
            csv_path = joinpath(out_table_dir, @sprintf("sunny_kpm_1d_%s_%gT.csv", qtag, B))
            sv_write_xy_csv(csv_path, "energy_meV,I_total,I_dispersive,I_flat,second_kernel_weight", disp.energy_meV, total, disp.intensity, flat.intensity, fill(r2, length(total)))
            push!(all_rows, (; B_T=B, qtag, csv_path))
        end
    end

    # Summary plot from the saved CSVs.
    fig = Figure(size=(1200, 750))
    for (iq, qtag) in enumerate(qtags)
        ax = Axis(fig[iq,1], xlabel="Energy (meV)", ylabel="arb.", title="Sunny KPM 1D $(qtag)")
        for B in fields
            pathcsv = joinpath(out_table_dir, @sprintf("sunny_kpm_1d_%s_%gT.csv", qtag, B))
            dat = readdlm(pathcsv, ',', skipstart=1)
            lines!(ax, dat[:,1], dat[:,2], label=@sprintf("%g T total", B))
        end
        axislegend(ax, position=:rt)
    end
    fig_path = joinpath(out_fig_dir, "sunny_kpm_1d_summary.png")
    save(fig_path, fig)
    @info "Saved KPM Sunny validation" fig_path
    return (; rows=all_rows, fig_path)
end

end # module
