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

function sv_system_size_controls(controls::Dict, section::AbstractString)
    sec = controls[section]

    if haskey(sec, "dims")
        dims = sv_dims_tuple(sec["dims"])
    elseif haskey(sec, "system_size")
        dims = sv_dims_tuple(sec["system_size"])
    else
        error("[$section] must define either dims or system_size")
    end

    if haskey(sec, "system_size")
        system_size = sv_dims_tuple(sec["system_size"])
        for i in 1:3
            system_size[i] >= dims[i] || error("[$section] system_size must be >= dims component-wise; got system_size=$system_size dims=$dims")
            system_size[i] % dims[i] == 0 || error("[$section] system_size must be an integer multiple of dims; got system_size=$system_size dims=$dims")
        end
        repeat_factor = (system_size[1] ÷ dims[1], system_size[2] ÷ dims[2], system_size[3] ÷ dims[3])

        if haskey(sec, "repeat_factor")
            configured_repeat = sv_dims_tuple(sec["repeat_factor"])
            implied_system_size = (dims[1] * configured_repeat[1], dims[2] * configured_repeat[2], dims[3] * configured_repeat[3])
            if implied_system_size != system_size
                @warn "Ignoring [$section].repeat_factor because [$section].system_size is present and inconsistent" section dims system_size configured_repeat implied_system_size repeat_factor
            end
        end
    else
        repeat_factor = haskey(sec, "repeat_factor") ? sv_dims_tuple(sec["repeat_factor"]) : (1, 1, 1)
        system_size = (dims[1] * repeat_factor[1], dims[2] * repeat_factor[2], dims[3] * repeat_factor[3])
    end

    return (; dims, repeat_factor, system_size)
end

function sv_sunny_transverse_gxy(controls::Dict)
    # Sunny's ssf_perp measure uses the magnetic moment tensor.  For H || c,
    # inelastic neutron intensity comes from transverse fluctuations, so setting
    # gxx = gyy = 0 can make the KPM spectrum vanish numerically.
    #
    # This is a Sunny neutron-intensity gauge for validation, not a fitted
    # physical gperp.  The effective flat/dispersive transverse intensity ratio
    # remains handled externally through gperp_ratio in sv_flat_neutron_weight.
    common = get(controls, "common", Dict{String,Any}())
    return Float64(get(common, "sunny_transverse_gxy", 1.0))
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
            "flat_to_dispersive_fraction",
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

function sv_csv_cell(x)
    if x isa Real
        return isfinite(Float64(x)) ? @sprintf("%.10g", Float64(x)) : string(x)
    elseif x === missing || x === nothing
        return ""
    else
        s = string(x)
        # Quote strings only when CSV syntax requires it.
        if findfirst(==(','), s) !== nothing || findfirst(==(Char(34)), s) !== nothing || findfirst(==(Char(10)), s) !== nothing || findfirst(==(Char(13)), s) !== nothing
            return "\"" * replace(s, "\"" => "\"\"") * "\""
        else
            return s
        end
    end
end

function sv_write_xy_csv(path::AbstractString, header::AbstractString, cols...)
    mkpath(dirname(path))
    isempty(cols) && error("sv_write_xy_csv requires at least one column")
    n = length(cols[1])
    for (j, c) in enumerate(cols)
        length(c) == n || error("CSV column $j has length $(length(c)); expected $n")
    end
    open(path, "w") do io
        println(io, header)
        for i in 1:n
            for (j, c) in enumerate(cols)
                j > 1 && print(io, ',')
                print(io, sv_csv_cell(c[i]))
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

function sv_set_gzz_all_sites!(sys, gzz::Real; gxy::Real=1.0)
    for site in eachsite(sys)
        G = sys.gs[site]
        Gm = Matrix(G)
        Gm .= 0.0
        Gm[1,1] = Float64(gxy)
        Gm[2,2] = Float64(gxy)
        Gm[3,3] = Float64(gzz)
        sys.gs[site] = SMatrix{3,3,Float64,9}(Gm)
    end
    return sys
end

function sv_build_effective_sunny_system(params, controls::Dict; component::Symbol=:dispersive, dims=(4,4,1), field_T::Real=0.0)
    units = sv_units()
    cryst = sv_effective_triangle_crystal(controls)
    g0 = component == :flat ? params.gzz2 : params.gzz
    gxy = sv_sunny_transverse_gxy(controls)
    # Start from a benign moment tensor and then overwrite every site below.
    # Keeping gxx = gyy nonzero is essential for ssf_perp neutron intensity.
    sys = System(cryst, [1 => Moment(s=Float64(controls["common"]["spin_S"]), g=1.0)], :dipole; dims=dims, seed=Int(controls["common"]["seed"]))
    sv_set_gzz_all_sites!(sys, g0; gxy=gxy)

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
    sizectl = sv_system_size_controls(controls, "largecell")
    dims = sizectl.dims
    repeat_factor = sizectl.repeat_factor
    @info "Sunny finite-size convention: large-cell" dims repeat_factor system_size=sizectl.system_size
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


# -----------------------------------------------------------------------------
# Experimental neutron 1D cuts and histogramming helpers
# -----------------------------------------------------------------------------

struct SVNeutronCut1D
    path::String
    Ei_meV::Float64
    temperature_K::Float64
    field_T::Float64
    qtag::String
    energy_meV::Vector{Float64}
    intensity::Vector{Float64}
    error::Vector{Float64}
    raw_intensity::Vector{Float64}
    background::Vector{Float64}
    background_level::Float64
    data_mode::Symbol
end

function sv_parse_filename_number_token(tok::AbstractString)
    return parse(Float64, replace(tok, "p" => "."))
end

function sv_parse_neutron_1d_filename(path::AbstractString)
    fname = basename(path)
    m = match(r"^yzgo_(\d+(?:p\d+)?)meV_(\d+(?:p\d+)?)K_(\d+(?:p\d+)?)T_Escan_(.+)_SYM\.dat$", fname)
    m === nothing && error("Filename does not match expected 1D scan pattern: $fname")
    return (;
        Ei_meV = sv_parse_filename_number_token(m.captures[1]),
        temperature_K = sv_parse_filename_number_token(m.captures[2]),
        field_T = sv_parse_filename_number_token(m.captures[3]),
        qtag = String(m.captures[4]),
    )
end

function sv_read_numeric_dat_matrix(path::AbstractString)
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            parts = split(s)
            vals = Float64[]
            good = true
            for part in parts
                v = tryparse(Float64, part)
                if v === nothing
                    good = false
                    break
                end
                push!(vals, v)
            end
            if good && !isempty(vals)
                push!(rows, vals)
            end
        end
    end
    isempty(rows) && error("No numeric rows found in neutron data file: $path")
    ncols = minimum(length.(rows))
    mat = zeros(Float64, length(rows), ncols)
    for (i, row) in enumerate(rows)
        mat[i, :] .= row[1:ncols]
    end
    return mat
end


# -----------------------------------------------------------------------------
# Analytical-cofit-style neutron background subtraction
# -----------------------------------------------------------------------------

struct SVRawNeutronScan1D
    path::String
    Ei_meV::Float64
    temperature_K::Float64
    field_T::Float64
    qtag::String
    energy_meV::Vector{Float64}
    intensity::Vector{Float64}
    error::Vector{Float64}
    H::Vector{Float64}
    K::Vector{Float64}
    L::Vector{Float64}
end

function sv_load_neutron_raw_scan_1d(path::AbstractString, controls::Dict)
    meta = sv_parse_neutron_1d_filename(path)
    mat = sv_read_numeric_dat_matrix(path)
    kc = controls["kpm"]

    # Match the analytical co-fit loader exactly.
    #
    # CNCS 1D scan file columns:
    #   1 Intensity
    #   2 Error
    #   3 DeltaE
    #   4 [0,K,0]
    #   5 [0,0,L]
    #   6 [H,0,0]
    ncols = size(mat, 2)
    ncols >= 6 || error("Expected at least 6 columns in $(basename(path)); got $ncols")

    intensity = Float64.(mat[:, 1])
    error     = Float64.(mat[:, 2])
    energy    = Float64.(mat[:, 3])
    K         = Float64.(mat[:, 4])
    L         = Float64.(mat[:, 5])
    H         = Float64.(mat[:, 6])

    idx = sortperm(energy)
    return SVRawNeutronScan1D(
        String(path),
        meta.Ei_meV,
        meta.temperature_K,
        meta.field_T,
        meta.qtag,
        energy[idx],
        intensity[idx],
        error[idx],
        H[idx],
        K[idx],
        L[idx],
    )
end

function sv_energy_window_mask(E::AbstractVector{<:Real}, windows::Vector{Tuple{Float64,Float64}})
    mask = falses(length(E))
    for (lo, hi) in windows
        mask .|= (E .>= lo) .& (E .<= hi)
    end
    return mask
end

function sv_background_fields_from_controls(controls::Dict)
    kc = controls["kpm"]
    vals = get(kc, "background_fields_T", Any[0.0, 9.0, 14.0])
    return Float64.(vals)
end

function sv_structured_residual_windows_from_controls(controls::Dict)
    kc = controls["kpm"]
    # Keep defaults identical to the analytical co-fit background model.
    d = Dict{String,Tuple{Float64,Float64}}(
        "0p33_0p33_0" => (1.675, 2.375),
        "0p5_0_0" => (1.825, 2.425),
    )
    if haskey(kc, "structured_residual_windows")
        d = Dict{String,Tuple{Float64,Float64}}()
        for (qtag, win) in kc["structured_residual_windows"]
            length(win) == 2 || error("structured_residual_windows.$qtag must have length 2")
            d[String(qtag)] = (Float64(win[1]), Float64(win[2]))
        end
    end
    return d
end

function sv_common_energy_grid(byfield::Dict{Float64,SVRawNeutronScan1D}, fields::Vector{Float64}; atol=1e-10)
    missing = [B for B in fields if !haskey(byfield, B)]
    isempty(missing) || error("Missing required background field scans: $missing")
    Eref = byfield[fields[1]].energy_meV
    for B in fields[2:end]
        E = byfield[B].energy_meV
        length(E) == length(Eref) || error("Energy grid length mismatch for q=$(byfield[B].qtag), B=$B")
        maximum(abs.(E .- Eref)) <= atol || error("Energy grids are not identical for q=$(byfield[B].qtag), B=$B")
    end
    return copy(Eref)
end

function sv_min_over_fields_background_raw(byfield::Dict{Float64,SVRawNeutronScan1D}, fields::Vector{Float64}; low_window=(0.0,0.75), high_threshold=2.5)
    E = sv_common_energy_grid(byfield, fields)
    I_by_field = [byfield[B].intensity for B in fields]
    Imin = similar(E)
    for i in eachindex(E)
        Imin[i] = minimum(I[i] for I in I_by_field)
    end
    lo, hi = low_window
    mask = ((E .>= lo) .& (E .<= hi)) .| (E .> high_threshold)
    return E[mask], Imin[mask]
end

function sv_sort_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || error("x and y must have same length")
    p = sortperm(x)
    return Float64.(x[p]), Float64.(y[p])
end

function sv_gaussian_smooth_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, sigma::Real)
    if sigma <= 0
        return Float64.(y)
    end
    ys = similar(Float64.(y))
    σ2 = Float64(sigma)^2
    for i in eachindex(x)
        wsum = 0.0
        ysum = 0.0
        xi = Float64(x[i])
        for j in eachindex(x)
            dx = Float64(x[j]) - xi
            w = exp(-0.5 * dx^2 / σ2)
            wsum += w
            ysum += w * Float64(y[j])
        end
        ys[i] = ysum / wsum
    end
    return ys
end

struct SVPchipInterpolator
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

function sv_pchip_endpoint_slope(h1::Float64, h2::Float64, d1::Float64, d2::Float64)
    m = ((2.0*h1 + h2)*d1 - h1*d2) / (h1 + h2)
    if sign(m) != sign(d1)
        return 0.0
    elseif sign(d1) != sign(d2) && abs(m) > abs(3.0*d1)
        return 3.0*d1
    else
        return m
    end
end

function sv_pchip_interpolator(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real})
    x, y = sv_sort_xy(xin, yin)
    n = length(x)
    n >= 2 || error("Need at least two points for PCHIP interpolation")
    any(diff(x) .<= 0) && error("PCHIP x-values must be strictly increasing")
    if n == 2
        d = (y[2] - y[1]) / (x[2] - x[1])
        return SVPchipInterpolator(x, y, [d, d])
    end
    h = diff(x)
    d = diff(y) ./ h
    m = zeros(Float64, n)
    m[1] = sv_pchip_endpoint_slope(h[1], h[2], d[1], d[2])
    m[n] = sv_pchip_endpoint_slope(h[end], h[end-1], d[end], d[end-1])
    for k in 2:n-1
        if d[k-1] == 0.0 || d[k] == 0.0 || sign(d[k-1]) != sign(d[k])
            m[k] = 0.0
        else
            w1 = 2.0*h[k] + h[k-1]
            w2 = h[k] + 2.0*h[k-1]
            m[k] = (w1 + w2) / (w1/d[k-1] + w2/d[k])
        end
    end
    return SVPchipInterpolator(x, y, m)
end

function (p::SVPchipInterpolator)(x0::Real)
    x = p.x; y = p.y; m = p.m
    n = length(x)
    j = searchsortedlast(x, Float64(x0))
    j = clamp(j, 1, n - 1)
    h = x[j+1] - x[j]
    t = (Float64(x0) - x[j]) / h
    h00 = 2.0*t^3 - 3.0*t^2 + 1.0
    h10 = t^3 - 2.0*t^2 + t
    h01 = -2.0*t^3 + 3.0*t^2
    h11 = t^3 - t^2
    return h00*y[j] + h10*h*m[j] + h01*y[j+1] + h11*h*m[j+1]
end

(p::SVPchipInterpolator)(xv::AbstractVector{<:Real}) = [p(x) for x in xv]

function sv_make_interpolated_background(Egrid::AbstractVector{<:Real}, Eraw::AbstractVector{<:Real}, Iraw::AbstractVector{<:Real}; smooth_sigma_meV::Real=0.0, interpolation_kind::Symbol=:pchip)
    Es, Is = sv_sort_xy(Eraw, Iraw)
    Is_smooth = sv_gaussian_smooth_xy(Es, Is, smooth_sigma_meV)
    bg = if interpolation_kind == :linear
        sv_interp1(Es, Is_smooth, Float64.(Egrid))
    elseif interpolation_kind == :pchip
        sv_pchip_interpolator(Es, Is_smooth)(Egrid)
    else
        error("Unknown interpolation_kind=$interpolation_kind. Use :pchip or :linear for Sunny validation background.")
    end
    return bg, Es, Is_smooth
end

function sv_ridge_linear_least_squares(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, err=nothing; ridge_lambda::Real=0.0, unpenalized_columns::Vector{Int}=Int[])
    Xf = Float64.(X); yf = Float64.(y)
    if err !== nothing
        σ = max.(Float64.(err), eps(Float64))
        sw = 1.0 ./ σ
        Xf = Xf .* sw
        yf = yf .* sw
    end
    A = transpose(Xf) * Xf
    b = transpose(Xf) * yf
    if ridge_lambda > 0
        penalty = Matrix{Float64}(I, size(A,1), size(A,2))
        for j in unpenalized_columns
            penalty[j,j] = 0.0
        end
        A .+= Float64(ridge_lambda) .* penalty
    end
    return A \ b
end

function sv_continuum_design_matrix(E::AbstractVector{<:Real}; model::Symbol=:power_tail, power::Real=3.0, energy_offset::Real=0.15, include_linear_tilt::Bool=false)
    Ef = Float64.(E)
    cols = Vector{Vector{Float64}}()
    push!(cols, ones(length(Ef)))
    if model == :power_tail
        push!(cols, (Ef .+ Float64(energy_offset)) .^ (-Float64(power)))
    elseif model == :exp_tail
        push!(cols, exp.(-Ef ./ Float64(power)))
    elseif model == :line
        push!(cols, Ef .- mean(Ef))
    else
        error("Unknown continuum model=$model")
    end
    if include_linear_tilt && model != :line
        push!(cols, Ef .- mean(Ef))
    end
    return hcat(cols...)
end

function sv_eval_continuum_model(E::AbstractVector{<:Real}, coeff::AbstractVector{<:Real}; model::Symbol=:power_tail, power::Real=3.0, energy_offset::Real=0.15, include_linear_tilt::Bool=false, center_energy::Real=0.0)
    Ef = Float64.(E)
    y = fill(Float64(coeff[1]), length(Ef))
    j = 2
    if model == :power_tail
        y .+= Float64(coeff[j]) .* (Ef .+ Float64(energy_offset)) .^ (-Float64(power)); j += 1
    elseif model == :exp_tail
        y .+= Float64(coeff[j]) .* exp.(-Ef ./ Float64(power)); j += 1
    elseif model == :line
        y .+= Float64(coeff[j]) .* (Ef .- Float64(center_energy)); j += 1
    else
        error("Unknown continuum model=$model")
    end
    if include_linear_tilt && model != :line
        y .+= Float64(coeff[j]) .* (Ef .- Float64(center_energy))
    end
    return y
end

function sv_fit_zeroT_continuum_baseline(Efit::AbstractVector{<:Real}, Ifit::AbstractVector{<:Real}; errfit=nothing, model::Symbol=:power_tail, power_grid=collect(0.5:0.05:8.0), exp_tau_grid=collect(0.25:0.025:3.0), energy_offset::Real=0.15, ridge_lambda::Real=0.0, include_linear_tilt::Bool=false, positive_tail::Bool=true)
    Ef = Float64.(Efit); If = Float64.(Ifit)
    center_energy = mean(Ef)
    scan_grid = model == :power_tail ? collect(power_grid) : model == :exp_tail ? collect(exp_tau_grid) : [NaN]
    best_score = Inf; best_coeff = Float64[]; best_param = NaN
    for param in scan_grid
        X = sv_continuum_design_matrix(Ef; model, power=isnan(param) ? 3.0 : param, energy_offset, include_linear_tilt)
        coeff = sv_ridge_linear_least_squares(X, If, errfit; ridge_lambda, unpenalized_columns=[1])
        if positive_tail && (model == :power_tail || model == :exp_tail) && coeff[2] < 0
            continue
        end
        pred = X * coeff
        resid = If .- pred
        score = if errfit !== nothing
            σ = max.(Float64.(errfit), eps(Float64))
            mean((resid ./ σ).^2)
        else
            mean(resid.^2)
        end
        if score < best_score
            best_score = score
            best_coeff = Float64.(coeff)
            best_param = isnan(param) ? 3.0 : Float64(param)
        end
    end
    if isempty(best_coeff)
        return sv_fit_zeroT_continuum_baseline(Efit, Ifit; errfit, model, power_grid, exp_tau_grid, energy_offset, ridge_lambda, include_linear_tilt, positive_tail=false)
    end
    baseline(Enew) = sv_eval_continuum_model(Enew, best_coeff; model, power=best_param, energy_offset, include_linear_tilt, center_energy)
    return baseline, best_coeff, best_param, best_score
end

function sv_structured_residual_points(byfield::Dict{Float64,SVRawNeutronScan1D}, residual_window::Tuple{Float64,Float64}; fit_window=(1.0,3.0), model::Symbol=:power_tail, power_grid=collect(0.5:0.05:8.0), exp_tau_grid=collect(0.25:0.025:3.0), energy_offset::Real=0.15, ridge_lambda::Real=0.0, include_linear_tilt::Bool=false, positive_tail::Bool=true, clip_negative_residuals::Bool=false)
    s0 = byfield[0.0]
    E = s0.energy_meV
    I0 = s0.intensity
    fit_lo, fit_hi = fit_window
    res_lo, res_hi = residual_window
    peak_mask = (E .>= res_lo) .& (E .<= res_hi)
    fit_mask = (E .>= fit_lo) .& (E .<= fit_hi) .& .!peak_mask
    any(fit_mask) || error("No points for 0 T continuum fit outside residual window")
    any(peak_mask) || error("No points inside residual window $residual_window")
    Efit = E[fit_mask]; Ifit = I0[fit_mask]; errfit = s0.error[fit_mask]
    baseline_fun, coeff, param, score = sv_fit_zeroT_continuum_baseline(Efit, Ifit; errfit, model, power_grid, exp_tau_grid, energy_offset, ridge_lambda, include_linear_tilt, positive_tail)
    Eres = E[peak_mask]
    baseline = baseline_fun(Eres)
    residual = I0[peak_mask] .- baseline
    if clip_negative_residuals
        residual = max.(residual, 0.0)
    end
    @printf("%s 0T continuum baseline: model=%s, param=%.4g, score=%.4g, coeff=%s\n", s0.qtag, String(model), param, score, repr(coeff))
    return Eres, residual, Efit, Ifit, baseline
end

function sv_make_analytical_background_model(qtag::String, byfield::Dict{Float64,SVRawNeutronScan1D}, controls::Dict)
    kc = controls["kpm"]
    fields = sv_background_fields_from_controls(controls)
    low_window = Tuple(Float64.(get(kc, "min_bg_low_window_meV", Any[0.0, 0.75])))
    high_threshold = Float64(get(kc, "min_bg_high_threshold_meV", 2.5))
    fit_window = Tuple(Float64.(get(kc, "structured_fit_window_meV", Any[1.0, 3.0])))
    final_smooth_sigma_meV = Float64(get(kc, "background_smooth_sigma_meV", 0.0))
    final_interp_kind = Symbol(get(kc, "background_interp_kind", "pchip"))
    residual_windows = sv_structured_residual_windows_from_controls(controls)
    zeroT_baseline_model = Symbol(get(kc, "zeroT_baseline_model", "power_tail"))
    energy_offset = Float64(get(kc, "zeroT_baseline_energy_offset_meV", 0.15))
    ridge_lambda = Float64(get(kc, "zeroT_baseline_ridge_lambda", 0.0))
    include_linear_tilt = Bool(get(kc, "zeroT_baseline_include_linear_tilt", false))
    positive_tail = Bool(get(kc, "zeroT_baseline_positive_tail", true))
    clip_negative_residuals = Bool(get(kc, "clip_negative_structured_residuals", false))

    Egrid = sv_common_energy_grid(byfield, fields)
    E_lowhigh, I_lowhigh = sv_min_over_fields_background_raw(byfield, fields; low_window, high_threshold)
    Eraw = copy(E_lowhigh); Iraw = copy(I_lowhigh)

    if haskey(residual_windows, qtag)
        Eres, residual, _, _, _ = sv_structured_residual_points(byfield, residual_windows[qtag]; fit_window, model=zeroT_baseline_model, energy_offset, ridge_lambda, include_linear_tilt, positive_tail, clip_negative_residuals)
        bg_without_residual, _, _ = sv_make_interpolated_background(Egrid, E_lowhigh, I_lowhigh; smooth_sigma_meV=final_smooth_sigma_meV, interpolation_kind=final_interp_kind)
        residual_base = sv_interp1(Float64.(Egrid), Float64.(bg_without_residual), Eres)
        residual_abs_bg = residual_base .+ residual
        append!(Eraw, Eres)
        append!(Iraw, residual_abs_bg)
    end

    bg, _, _ = sv_make_interpolated_background(Egrid, Eraw, Iraw; smooth_sigma_meV=final_smooth_sigma_meV, interpolation_kind=final_interp_kind)
    return Float64.(bg)
end

function sv_constant_tail_background(E::AbstractVector, I::AbstractVector, controls::Dict)
    kc = controls["kpm"]
    windows = get(kc, "tail_background_windows_meV", Any[[0.0, 0.75], [2.5, 4.0]])
    vals = Float64[]
    for win in windows
        length(win) == 2 || error("Each tail_background_windows_meV entry must have length 2")
        lo = Float64(win[1]); hi = Float64(win[2])
        for (e, y) in zip(E, I)
            if isfinite(e) && isfinite(y) && lo <= e <= hi
                push!(vals, Float64(y))
            end
        end
    end
    isempty(vals) && return zeros(Float64, length(E))
    method = Symbol(get(kc, "tail_background_statistic", "median"))
    bg = method == :mean ? mean(vals) : method == :median ? median(vals) : error("Unknown tail_background_statistic=$method")
    return fill(bg, length(E))
end

function sv_make_corrected_cut(raw::SVRawNeutronScan1D, bg::AbstractVector, controls::Dict)
    mode = Symbol(get(controls["kpm"], "data_mode", "tail_bgsub"))
    length(bg) == length(raw.energy_meV) || error("Background length mismatch for $(raw.path)")
    I = if mode == :raw || mode == :file_intensity
        Float64.(raw.intensity)
    elseif mode in (:tail_bgsub, :spline_bgsub, :bgsub, :analytical_bgsub)
        Float64.(raw.intensity) .- Float64.(bg)
    else
        error("Unknown kpm data_mode=$mode")
    end
    return SVNeutronCut1D(
        raw.path,
        raw.Ei_meV,
        raw.temperature_K,
        raw.field_T,
        raw.qtag,
        copy(raw.energy_meV),
        I,
        copy(raw.error),
        copy(raw.intensity),
        mode == :raw || mode == :file_intensity ? zeros(Float64, length(raw.energy_meV)) : Float64.(bg),
        mode == :raw || mode == :file_intensity ? 0.0 : mean(Float64.(bg)),
        mode,
    )
end

function sv_make_corrected_cuts_for_qtag(qtag::String, byfield::Dict{Float64,SVRawNeutronScan1D}, controls::Dict)
    mode = Symbol(get(controls["kpm"], "data_mode", "tail_bgsub"))
    bg_by_field = Dict{Float64,Vector{Float64}}()
    if mode == :raw || mode == :file_intensity
        for (B, raw) in byfield
            bg_by_field[B] = zeros(Float64, length(raw.energy_meV))
        end
    elseif mode in (:tail_bgsub, :spline_bgsub, :bgsub, :analytical_bgsub)
        bg = sv_make_analytical_background_model(qtag, byfield, controls)
        for (B, raw) in byfield
            bg_by_field[B] = copy(bg)
        end
    elseif mode == :constant_tail_bgsub
        for (B, raw) in byfield
            bg_by_field[B] = sv_constant_tail_background(raw.energy_meV, raw.intensity, controls)
        end
    else
        error("Unknown kpm data_mode=$mode")
    end
    return Dict(B => sv_make_corrected_cut(raw, bg_by_field[B], controls) for (B, raw) in byfield)
end

function sv_load_kpm_experimental_cuts(repo_root::AbstractString, controls::Dict)
    kc = controls["kpm"]
    dir = sv_repo_path(repo_root, controls["paths"]["neutron_1d_dir"])
    isdir(dir) || error("Could not find neutron 1D directory: $dir")

    target_Ei = Float64(get(kc, "Ei_meV", 4.65))
    target_T = Float64(get(kc, "temperature_K", 0.07))
    target_fields = Float64.(controls["common"]["fields_T"])
    background_fields = sv_background_fields_from_controls(controls)
    needed_fields = unique(vcat(target_fields, background_fields))
    qtags = Set(String.(kc["qtags"]))

    raw_by_q = Dict{String,Dict{Float64,SVRawNeutronScan1D}}()
    for fname in sort(readdir(dir))
        endswith(lowercase(fname), ".dat") || continue
        path = joinpath(dir, fname)
        meta = try
            sv_parse_neutron_1d_filename(path)
        catch
            continue
        end
        sv_nearly_equal(meta.Ei_meV, target_Ei; atol=1e-3) || continue
        sv_nearly_equal(meta.temperature_K, target_T; atol=1e-3) || continue
        meta.qtag in qtags || continue
        any(B -> sv_nearly_equal(meta.field_T, B; atol=1e-3), needed_fields) || continue
        raw = sv_load_neutron_raw_scan_1d(path, controls)
        if !haskey(raw_by_q, raw.qtag)
            raw_by_q[raw.qtag] = Dict{Float64,SVRawNeutronScan1D}()
        end
        raw_by_q[raw.qtag][raw.field_T] = raw
    end

    cuts = SVNeutronCut1D[]
    for qtag in sort(collect(qtags))
        haskey(raw_by_q, qtag) || error("No raw scans found for qtag=$qtag")
        corrected = sv_make_corrected_cuts_for_qtag(qtag, raw_by_q[qtag], controls)
        for B in target_fields
            haskey(corrected, B) || error("No corrected target-field cut for qtag=$qtag, B=$B")
            push!(cuts, corrected[B])
        end
    end

    sort!(cuts; by = c -> (c.field_T, c.qtag))
    return cuts
end

function sv_nearly_equal(a::Real, b::Real; atol::Real=1e-6)
    return abs(Float64(a) - Float64(b)) <= Float64(atol)
end

function sv_find_cut(cuts::Vector{SVNeutronCut1D}, field_T::Real, qtag::AbstractString)
    for cut in cuts
        if sv_nearly_equal(cut.field_T, field_T; atol=1e-3) && cut.qtag == qtag
            return cut
        end
    end
    return nothing
end

function sv_energy_bin_edges(E::AbstractVector)
    n = length(E)
    n >= 2 || error("Need at least two energy points to infer bin edges")
    edges = zeros(Float64, n + 1)
    for i in 2:n
        edges[i] = 0.5 * (E[i-1] + E[i])
    end
    edges[1] = E[1] - 0.5 * (E[2] - E[1])
    edges[end] = E[end] + 0.5 * (E[end] - E[end-1])
    return edges
end

function sv_model_to_experimental_energy_grid(
    E_model::AbstractVector,
    I_model::AbstractVector,
    E_exp::AbstractVector;
    mode::Symbol = :bin_average,
)
    Em = Float64.(E_model)
    Im = Float64.(I_model)
    Ee = Float64.(E_exp)

    idx = sortperm(Em)
    Em = Em[idx]
    Im = Im[idx]

    if mode == :interpolate
        return sv_interp1(Em, Im, Ee)
    elseif mode == :bin_average
        edges = sv_energy_bin_edges(Ee)
        out = fill(NaN, length(Ee))
        for i in eachindex(Ee)
            lo = edges[i]
            hi = edges[i+1]
            mask = (Em .>= lo) .& (Em .< hi) .& isfinite.(Im)
            if any(mask)
                out[i] = mean(Im[mask])
            else
                out[i] = sv_interp1(Em, Im, [Ee[i]])[1]
            end
        end
        return out
    else
        error("Unknown histogram/interpolation mode = $mode")
    end
end

function sv_neutron_scale(params, controls::Dict)
    kc = controls["kpm"]
    if get(kc, "use_neutron_global_scale_from_best_fit", true)
        hasproperty(params, :neutron_global_scale) || error("Canonical params do not include neutron_global_scale")
        return Float64(params.neutron_global_scale)
    else
        return Float64(get(kc, "manual_neutron_global_scale", 1.0))
    end
end

function sv_flat_neutron_weight(params, controls::Dict)
    r2 = sv_second_kernel_weight(params, controls)
    kc = controls["kpm"]
    if get(kc, "include_gperp_ratio_intensity_weight", true)
        return r2 * params.gperp_ratio^2
    else
        return r2
    end
end

function sv_write_neutron_cut_inventory(path::AbstractString, cuts::Vector{SVNeutronCut1D})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "path,Ei_meV,temperature_K,field_T,qtag,n_E,E_min_meV,E_max_meV,I_min,I_max,rawI_min,rawI_max,bg_min,bg_max,data_mode,background_mean")
        for c in cuts
            println(io, join([
                c.path,
                c.Ei_meV,
                c.temperature_K,
                c.field_T,
                c.qtag,
                length(c.energy_meV),
                minimum(c.energy_meV),
                maximum(c.energy_meV),
                minimum(c.intensity),
                maximum(c.intensity),
                minimum(c.raw_intensity),
                maximum(c.raw_intensity),
                minimum(c.background),
                maximum(c.background),
                String(c.data_mode),
                c.background_level,
            ], ","))
        end
    end
end

function sv_check_neutron_cut_loading(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    cuts = sv_load_kpm_experimental_cuts(repo_root, controls)
    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    mkpath(out_table_dir)
    inventory_path = joinpath(out_table_dir, "neutron_1d_cut_inventory.csv")
    sv_write_neutron_cut_inventory(inventory_path, cuts)

    println("Loaded neutron 1D cuts for Sunny validation")
    println("------------------------------------------")
    for c in cuts
        @printf("B=%6.3f T  q=%-14s  n=%4d  E=[%7.3f,%7.3f] meV  I=[% .4g,% .4g]  mode=%s  bg=[% .4g,% .4g] mean=% .6g\n",
            c.field_T, c.qtag, length(c.energy_meV),
            minimum(c.energy_meV), maximum(c.energy_meV),
            minimum(c.intensity), maximum(c.intensity),
            String(c.data_mode), minimum(c.background), maximum(c.background), c.background_level)
    end
    println()
    println("Wrote inventory:")
    println(inventory_path)
    return (; cuts, inventory_path)
end


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
    sizectl = sv_system_size_controls(controls, "kpm")
    dims = sizectl.dims
    repeat_factor = sizectl.repeat_factor
    include_exchange = get(kc, "include_exchange_disorder", true)
    include_gzz = get(kc, "include_gzz_disorder", true)
    maxiters = Int(kc["maxiters"])

    base = sv_build_effective_sunny_system(params, controls; component, dims, field_T)
    sys = base.sys
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
    res = intensities(swt, qs; energies, kernel)
    raw = sv_try_extract_sunny_intensity(res)

    I = if raw isa AbstractVector
        Float64.(raw)
    elseif raw isa AbstractMatrix
        # Common Sunny layouts are (nq, nE) or (nE, nq).  Since we supply one
        # q-point, choose the orientation that matches the energy axis.
        r, c = size(raw)
        if r == length(energies)
            vec(mean(Float64.(raw); dims=2))
        elseif c == length(energies)
            vec(mean(Float64.(raw); dims=1))
        else
            vec(Float64.(raw))[1:min(length(raw), length(energies))]
        end
    else
        error("Unsupported Sunny intensity container type: $(typeof(raw))")
    end

    if length(I) != length(energies)
        tmp = fill(NaN, length(energies))
        n = min(length(I), length(energies))
        tmp[1:n] .= I[1:n]
        I = tmp
    end
    return (; energy_meV=energies, intensity=I, result=res)
end


function sv_kpm_neutron_scale_mode(controls::Dict)
    kc = controls["kpm"]
    if haskey(kc, "neutron_scale_mode")
        return Symbol(kc["neutron_scale_mode"])
    end

    # Backward-compatible interpretation of the older controls.
    if get(kc, "use_neutron_global_scale_from_best_fit", true)
        return :best_fit
    else
        return :manual
    end
end

function sv_kpm_neutron_scale_scope(controls::Dict)
    kc = controls["kpm"]
    return Symbol(get(kc, "neutron_scale_scope", "global"))
end

function sv_kpm_neutron_scale_mask(E::AbstractVector, y::AbstractVector, x::AbstractVector, err::AbstractVector, controls::Dict)
    kc = controls["kpm"]
    mask = isfinite.(E) .& isfinite.(y) .& isfinite.(x) .& isfinite.(err)

    if haskey(kc, "neutron_scale_fit_window_meV")
        win = Float64.(kc["neutron_scale_fit_window_meV"])
        length(win) == 2 || error("kpm.neutron_scale_fit_window_meV must have two entries")
        mask .&= (E .>= win[1]) .& (E .<= win[2])
    end

    if get(kc, "neutron_scale_positive_experiment_only", false)
        mask .&= y .> 0
    end

    if get(kc, "neutron_scale_positive_model_only", false)
        mask .&= x .> 0
    end

    return mask
end

function sv_kpm_scale_weights(err::AbstractVector, controls::Dict)
    kc = controls["kpm"]
    if get(kc, "neutron_scale_use_uncertainties", true)
        e = Float64.(err)
        finite_positive = e[isfinite.(e) .& (e .> 0)]
        floor = isempty(finite_positive) ? 1.0 : minimum(finite_positive)
        e = map(v -> (isfinite(v) && v > 0) ? v : floor, e)
        return 1.0 ./ (e .^ 2)
    else
        return ones(Float64, length(err))
    end
end

function sv_best_positive_scale_least_squares(y::AbstractVector, x::AbstractVector, err::AbstractVector, E::AbstractVector, controls::Dict; fallback::Float64=1.0)
    mask = sv_kpm_neutron_scale_mask(E, y, x, err, controls)
    if !any(mask)
        @warn "No finite points available for least-squares neutron scale; using fallback" fallback
        return fallback
    end
    w = sv_kpm_scale_weights(err, controls)
    xx = Float64.(x[mask])
    yy = Float64.(y[mask])
    ww = Float64.(w[mask])
    denom = sum(ww .* xx .* xx)
    min_model_power = Float64(get(controls["kpm"], "neutron_scale_min_model_power", 1e-30))
    if !(isfinite(denom) && denom > min_model_power)
        @warn "Degenerate or nearly-zero Sunny model for least-squares neutron scale; using fallback" fallback denom min_model_power maximum_abs_model=maximum(abs.(xx))
        return fallback
    end
    s = sum(ww .* xx .* yy) / denom
    if get(controls["kpm"], "neutron_scale_nonnegative", true)
        s = max(0.0, s)
    end
    max_abs_scale = Float64(get(controls["kpm"], "neutron_scale_max_abs", Inf))
    if !(isfinite(s))
        @warn "Non-finite least-squares neutron scale; using fallback" fallback scale=s
        return fallback
    elseif abs(s) > max_abs_scale
        @warn "Least-squares neutron scale exceeded configured guard; using fallback" fallback scale=s max_abs_scale denom maximum_abs_model=maximum(abs.(xx))
        return fallback
    end
    return Float64(s)
end

function sv_best_positive_scale_max_match(y::AbstractVector, x::AbstractVector, err::AbstractVector, E::AbstractVector, controls::Dict; fallback::Float64=1.0)
    mask = sv_kpm_neutron_scale_mask(E, y, x, err, controls)
    if !any(mask)
        @warn "No finite points available for max-match neutron scale; using fallback" fallback
        return fallback
    end
    yy = Float64.(y[mask])
    xx = Float64.(x[mask])
    ymax = maximum(yy)
    xmax = maximum(xx)
    min_model_peak = Float64(get(controls["kpm"], "neutron_scale_min_model_peak", 1e-15))
    if !(isfinite(ymax) && isfinite(xmax) && abs(xmax) > min_model_peak)
        @warn "Degenerate or nearly-zero model/data for max-match neutron scale; using fallback" fallback ymax xmax min_model_peak
        return fallback
    end
    s = ymax / xmax
    if get(controls["kpm"], "neutron_scale_nonnegative", true)
        s = max(0.0, s)
    end
    max_abs_scale = Float64(get(controls["kpm"], "neutron_scale_max_abs", Inf))
    if !(isfinite(s))
        @warn "Non-finite max-match neutron scale; using fallback" fallback scale=s
        return fallback
    elseif abs(s) > max_abs_scale
        @warn "Max-match neutron scale exceeded configured guard; using fallback" fallback scale=s max_abs_scale ymax xmax
        return fallback
    end
    return Float64(s)
end

function sv_compute_neutron_scale(mode::Symbol, y::AbstractVector, x::AbstractVector, err::AbstractVector, E::AbstractVector, controls::Dict; fallback::Float64=1.0)
    if mode == :best_fit
        return fallback
    elseif mode == :manual
        return fallback
    elseif mode == :least_squares
        return sv_best_positive_scale_least_squares(y, x, err, E, controls; fallback)
    elseif mode == :max_match
        return sv_best_positive_scale_max_match(y, x, err, E, controls; fallback)
    else
        error("Unknown kpm.neutron_scale_mode = $mode. Use best_fit, manual, least_squares, or max_match.")
    end
end

function sv_run_kpm_1d(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    print_canonical_model_parameters(params)
    @info "Sunny validation: KPM 1D spin-wave" sunny_version=sv_try_pkgversion(Sunny) params_path=path

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    out_fig_dir = sv_repo_path(repo_root, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir); mkpath(out_fig_dir)

    cuts = sv_load_kpm_experimental_cuts(repo_root, controls)
    inventory_path = joinpath(out_table_dir, "sunny_kpm_neutron_1d_cut_inventory.csv")
    sv_write_neutron_cut_inventory(inventory_path, cuts)

    fields = Float64.(controls["common"]["fields_T"])
    qtags = String.(controls["kpm"]["qtags"])
    fallback_neutron_scale = sv_neutron_scale(params, controls)
    scale_mode = sv_kpm_neutron_scale_mode(controls)
    scale_scope = sv_kpm_neutron_scale_scope(controls)
    flat_weight = sv_flat_neutron_weight(params, controls)
    hist_mode = Symbol(get(controls["kpm"], "histogram_mode", "bin_average"))
    sunny_transverse_gxy = sv_sunny_transverse_gxy(controls)
    kpm_sizectl = sv_system_size_controls(controls, "kpm")
    flat_to_dispersive_fraction = sv_second_kernel_weight(params, controls)

    @info "KPM comparison convention" scale_mode scale_scope fallback_neutron_scale flat_to_dispersive_fraction flat_weight histogram_mode=hist_mode sunny_transverse_gxy dims=kpm_sizectl.dims repeat_factor=kpm_sizectl.repeat_factor system_size=kpm_sizectl.system_size data_cuts=length(cuts)

    # First compute all unscaled Sunny spectra on the experimental energy grids.
    # The neutron intensity scale can then be either fixed, fitted globally, or
    # fitted per cut without rerunning Sunny/KPM.
    cut_results = NamedTuple[]
    for B in fields
        for qtag in qtags
            cut = sv_find_cut(cuts, B, qtag)
            cut === nothing && error("No experimental 1D cut found for B=$B T qtag=$qtag. Check kpm Ei/T/field/qtag controls.")

            @info "Computing Sunny KPM cut and histogramming to experiment" B_T=B qtag=qtag nE_exp=length(cut.energy_meV)
            disp = sv_kpm_component_spectrum(params, controls; component=:dispersive, field_T=B, qtag)
            flat = sv_kpm_component_spectrum(params, controls; component=:flat, field_T=B, qtag)

            Idisp_grid = sv_model_to_experimental_energy_grid(disp.energy_meV, disp.intensity, cut.energy_meV; mode=hist_mode)
            Iflat_grid = sv_model_to_experimental_energy_grid(flat.energy_meV, flat.intensity, cut.energy_meV; mode=hist_mode)
            Itotal_unscaled = Idisp_grid .+ flat_weight .* Iflat_grid

            push!(cut_results, (;
                B_T = B,
                qtag = qtag,
                cut = cut,
                Idisp_grid = Idisp_grid,
                Iflat_grid = Iflat_grid,
                Itotal_unscaled = Itotal_unscaled,
            ))
        end
    end

    global_scale = fallback_neutron_scale
    if scale_scope == :global
        Ecat = reduce(vcat, [r.cut.energy_meV for r in cut_results])
        ycat = reduce(vcat, [r.cut.intensity for r in cut_results])
        ecat = reduce(vcat, [r.cut.error for r in cut_results])
        xcat = reduce(vcat, [r.Itotal_unscaled for r in cut_results])
        global_scale = sv_compute_neutron_scale(scale_mode, ycat, xcat, ecat, Ecat, controls; fallback=fallback_neutron_scale)
    elseif scale_scope == :per_cut
        # Computed below for each cut.
    else
        error("Unknown kpm.neutron_scale_scope = $scale_scope. Use global or per_cut.")
    end

    @info "KPM neutron scale selected" scale_mode scale_scope global_scale fallback_neutron_scale

    all_rows = NamedTuple[]
    for r in cut_results
        cut = r.cut
        neutron_scale = if scale_scope == :per_cut
            sv_compute_neutron_scale(scale_mode, cut.intensity, r.Itotal_unscaled, cut.error, cut.energy_meV, controls; fallback=fallback_neutron_scale)
        else
            global_scale
        end

        Itotal_scaled = neutron_scale .* r.Itotal_unscaled
        Idisp_scaled = neutron_scale .* r.Idisp_grid
        Iflat_scaled = neutron_scale .* flat_weight .* r.Iflat_grid
        residual = cut.intensity .- Itotal_scaled

        csv_path = joinpath(out_table_dir, @sprintf("sunny_kpm_1d_%s_%gT_vs_exp.csv", r.qtag, r.B_T))
        sv_write_xy_csv(csv_path,
            "energy_meV,I_exp,Ierr_exp,I_total_scaled,I_disp_scaled,I_flat_scaled,I_total_unscaled,I_disp_unscaled,I_flat_unweighted,residual,raw_intensity_exp,background,neutron_scale,neutron_scale_mode,neutron_scale_scope,fallback_neutron_global_scale,flat_weight,flat_to_dispersive_fraction,gperp_ratio,sunny_dims_x,sunny_dims_y,sunny_dims_z,sunny_repeat_x,sunny_repeat_y,sunny_repeat_z,sunny_system_size_x,sunny_system_size_y,sunny_system_size_z",
            cut.energy_meV,
            cut.intensity,
            cut.error,
            Itotal_scaled,
            Idisp_scaled,
            Iflat_scaled,
            r.Itotal_unscaled,
            r.Idisp_grid,
            r.Iflat_grid,
            residual,
            cut.raw_intensity,
            cut.background,
            fill(neutron_scale, length(cut.energy_meV)),
            fill(String(scale_mode), length(cut.energy_meV)),
            fill(String(scale_scope), length(cut.energy_meV)),
            fill(fallback_neutron_scale, length(cut.energy_meV)),
            fill(flat_weight, length(cut.energy_meV)),
            fill(flat_to_dispersive_fraction, length(cut.energy_meV)),
            fill(params.gperp_ratio, length(cut.energy_meV)),
            fill(kpm_sizectl.dims[1], length(cut.energy_meV)),
            fill(kpm_sizectl.dims[2], length(cut.energy_meV)),
            fill(kpm_sizectl.dims[3], length(cut.energy_meV)),
            fill(kpm_sizectl.repeat_factor[1], length(cut.energy_meV)),
            fill(kpm_sizectl.repeat_factor[2], length(cut.energy_meV)),
            fill(kpm_sizectl.repeat_factor[3], length(cut.energy_meV)),
            fill(kpm_sizectl.system_size[1], length(cut.energy_meV)),
            fill(kpm_sizectl.system_size[2], length(cut.energy_meV)),
            fill(kpm_sizectl.system_size[3], length(cut.energy_meV)),
        )
        push!(all_rows, (; B_T=r.B_T, qtag=r.qtag, csv_path, neutron_scale))
    end

    fig = Figure(size=(1300, 850))
    ylims = get(controls["kpm"], "plot_ylim", nothing)
    for (iq, qtag) in enumerate(qtags)
        for (iB, B) in enumerate(fields)
            ax = Axis(fig[iq,iB], xlabel="Energy (meV)", ylabel="Intensity", title=@sprintf("%s, %g T", qtag, B))
            pathcsv = joinpath(out_table_dir, @sprintf("sunny_kpm_1d_%s_%gT_vs_exp.csv", qtag, B))
            dat = readdlm(pathcsv, ',', skipstart=1)
            scatter!(ax, dat[:,1], dat[:,2], markersize=6, label="experiment")
            lines!(ax, dat[:,1], dat[:,4], label="Sunny total")
            lines!(ax, dat[:,1], dat[:,5], linestyle=:dash, label="disp.")
            lines!(ax, dat[:,1], dat[:,6], linestyle=:dot, label="flat")
            if ylims !== nothing && length(ylims) == 2
                ylims!(ax, Float64(ylims[1]), Float64(ylims[2]))
            end
            iq == 1 && iB == length(fields) && axislegend(ax, position=:rt)
        end
    end
    fig_path = joinpath(out_fig_dir, "sunny_kpm_1d_vs_experiment.png")
    save(fig_path, fig)
    @info "Saved KPM Sunny validation" fig_path inventory_path scale_mode scale_scope global_scale fallback_neutron_scale
    return (; rows=all_rows, fig_path, inventory_path, cuts, scale_mode, scale_scope, global_scale, fallback_neutron_scale)
end

end # module
