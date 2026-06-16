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


# -----------------------------------------------------------------------------
# Yb3+ magnetic form factor and physical reciprocal-space helpers
# -----------------------------------------------------------------------------

function sv_radial_integral_j0(s::Real, A::Real, a::Real, B::Real, b::Real, C::Real, c::Real, D::Real)
    x = Float64(s)^2
    return A * exp(-a * x) + B * exp(-b * x) + C * exp(-c * x) + D
end

function sv_radial_integral_jL_nonzero(s::Real, A::Real, a::Real, B::Real, b::Real, C::Real, c::Real, D::Real)
    x = Float64(s)^2
    return (A * exp(-a * x) + B * exp(-b * x) + C * exp(-c * x) + D) * x
end

function sv_yb3_j0(Q_Ainv::Real)
    # Same coefficients/convention as the analytical 2D polarized-state model.
    s = Float64(Q_Ainv) / (4.0 * pi)
    return sv_radial_integral_j0(s, 0.0416, 16.0949, 0.2849, 7.8341, 0.6961, 2.6725, -0.0229)
end

function sv_yb3_j2(Q_Ainv::Real)
    # Same coefficients/convention as the analytical 2D polarized-state model.
    s = Float64(Q_Ainv) / (4.0 * pi)
    return sv_radial_integral_jL_nonzero(s, 0.1570, 18.5553, 0.8484, 6.5403, 0.8880, 2.0367, 0.0318)
end

function sv_yb3_form_factor(Q_Ainv::Real; include_j2::Bool=true, c2::Real=0.75)
    f = include_j2 ? sv_yb3_j0(Q_Ainv) + Float64(c2) * sv_yb3_j2(Q_Ainv) : sv_yb3_j0(Q_Ainv)
    return isfinite(f) ? f : 1.0
end

function sv_form_factor_controls(controls::Dict)
    ff = get(controls, "neutron_form_factor", Dict{String,Any}())
    ff isa Dict || (ff = Dict{String,Any}())

    enabled = Bool(get(ff, "enabled", false))

    # Prefer Sunny's built-in FormFactor table when computing Sunny intensities.
    # The manual Yb3+ implementation is retained only as a diagnostic/fallback
    # so we do not accidentally double-apply a magnetic form factor.
    source = Symbol(get(ff, "source", "sunny_builtin"))
    ion = String(get(ff, "ion", "Yb3"))
    candidate_ions = if haskey(ff, "candidate_ions")
        String.(ff["candidate_ions"])
    else
        unique(String[ion, replace(ion, "+"=>""), "Yb3", "Yb3+"])
    end

    manual_include_j2 = Bool(get(ff, "manual_include_j2", get(ff, "include_j2", true)))
    manual_c2 = Float64(get(ff, "manual_j2_coefficient", get(ff, "j2_coefficient", 0.75)))
    manual_apply_as = Symbol(get(ff, "manual_apply_as", get(ff, "apply_as", "intensity_squared")))

    on_error = Symbol(get(ff, "on_error", "error"))
    normalize_to = Symbol(get(ff, "normalize_to", "absolute_Q"))

    return (; enabled, source, ion, candidate_ions,
        include_j2=manual_include_j2, c2=manual_c2, apply_as=manual_apply_as,
        on_error, normalize_to)
end

function sv_builtin_formfactor_pairs(controls::Dict)
    ff = sv_form_factor_controls(controls)
    if !ff.enabled || ff.source != :sunny_builtin
        return nothing
    end

    last_err = nothing
    for ion in ff.candidate_ions
        try
            return [1 => FormFactor(ion)]
        catch err
            last_err = err
        end
    end

    msg = "Could not construct Sunny FormFactor for candidate ion labels $(ff.candidate_ions). " *
          "Try changing [neutron_form_factor].ion / candidate_ions in configs/sunny_validation_controls.toml."
    if ff.on_error in (:none, :off, :disable, :warn_disable)
        @warn msg exception=last_err
        return nothing
    else
        error(msg * " Last error: $(last_err)")
    end
end

function sv_sunny_measure(sys, controls::Dict)
    # The preferred path is Sunny-native form factors inside ssf_perp. This lets
    # Sunny handle the neutron measurement convention and prevents manual
    # |f(Q)|^2 double counting.
    pairs = sv_builtin_formfactor_pairs(controls)
    if pairs === nothing
        return ssf_perp(sys)
    else
        return ssf_perp(sys; formfactors=pairs)
    end
end

function sv_form_factor_lattice_controls(controls::Dict)
    ff = get(controls, "neutron_form_factor", Dict{String,Any}())
    lat = ff isa Dict ? get(ff, "lattice", Dict{String,Any}()) : Dict{String,Any}()
    lat isa Dict || (lat = Dict{String,Any}())
    common = get(controls, "common", Dict{String,Any}())
    a = Float64(get(lat, "a_A", get(common, "lattice_a_angstrom", 3.376)))
    c = Float64(get(lat, "c_A", get(common, "lattice_c_angstrom", 21.96)))
    gamma_deg = Float64(get(lat, "gamma_deg", 120.0))
    if abs(gamma_deg - 120.0) > 1e-6
        @warn "Sunny validation form-factor Q conversion currently assumes the analytical hexagonal/triangular gamma=120 convention" gamma_deg
    end
    return (; a_A=a, c_A=c, gamma_deg)
end

function sv_rlu_basis_matrix_for_form_factor(controls::Dict)
    lat = sv_form_factor_lattice_controls(controls)
    astar = 4.0 * pi / (sqrt(3.0) * lat.a_A)
    cstar = 2.0 * pi / lat.c_A

    # Same reciprocal basis as the analytical model:
    # Q = H a* + K b* + L c*, with a* and b* 120-degree hex/triangular duals.
    b1 = [ astar, 0.0, 0.0 ]
    b2 = [ -0.5 * astar, 0.5 * sqrt(3.0) * astar, 0.0 ]
    b3 = [ 0.0, 0.0, cstar ]
    return hcat(b1, b2, b3)
end

function sv_qcart_Ainv(q_hkl::AbstractVector{<:Real}, controls::Dict)
    Bmat = sv_rlu_basis_matrix_for_form_factor(controls)
    return Bmat * collect(Float64, q_hkl)
end

function sv_qmag_Ainv(q_hkl::AbstractVector{<:Real}, controls::Dict)
    return norm(sv_qcart_Ainv(q_hkl, controls))
end

function sv_form_factor_amplitude(q_hkl::AbstractVector{<:Real}, controls::Dict)
    ff = sv_form_factor_controls(controls)
    if !ff.enabled
        return 1.0
    elseif ff.source == :sunny_builtin
        # Sunny applies the form factor internally through ssf_perp(; formfactors).
        # Return unity here so diagnostics/CSV columns do not imply an additional
        # post-processing weight. Use source="manual_yb3" to print/apply this
        # analytical fallback explicitly.
        return 1.0
    elseif ff.source != :manual_yb3
        @warn "Unsupported neutron_form_factor source; using no manual post-factor" source=ff.source
        return 1.0
    end

    ion_norm = lowercase(replace(ff.ion, " "=>""))
    if !(ion_norm in ("yb3+", "yb3", "ybiii"))
        @warn "Unsupported manual neutron_form_factor ion; using no form factor" ion=ff.ion
        return 1.0
    end
    Q = sv_qmag_Ainv(q_hkl, controls)
    return sv_yb3_form_factor(Q; include_j2=ff.include_j2, c2=ff.c2)
end

function sv_form_factor_intensity_weight(q_hkl::AbstractVector{<:Real}, controls::Dict)
    ff = sv_form_factor_controls(controls)
    (!ff.enabled || ff.source == :sunny_builtin) && return 1.0
    f = sv_form_factor_amplitude(q_hkl, controls)
    if ff.apply_as == :intensity_squared
        return f^2
    elseif ff.apply_as == :amplitude
        return f
    elseif ff.apply_as == :none
        return 1.0
    else
        error("Unsupported [neutron_form_factor].manual_apply_as=$(ff.apply_as). Use intensity_squared, amplitude, or none.")
    end
end

function sv_apply_form_factor_to_intensity(I::AbstractMatrix, qs::Vector{Vector{Float64}}, controls::Dict)
    ff = sv_form_factor_controls(controls)
    size(I, 2) == length(qs) || error("Form factor application expected $(length(qs)) q columns but intensity has $(size(I,2))")

    # Always compute |Q| for diagnostics/CSV output.  For source = sunny_builtin,
    # Sunny applies the magnetic form factor inside ssf_perp(; formfactors), so
    # there is no additional post-processing weight here.  For source =
    # manual_yb3, apply the manual factor as a fallback.
    qmag = [sv_qmag_Ainv(q, controls) for q in qs]

    if !ff.enabled || ff.source == :sunny_builtin
        return Matrix{Float64}(I), ones(Float64, length(qs)), ones(Float64, length(qs)), qmag
    end

    amps = [sv_form_factor_amplitude(q, controls) for q in qs]
    weights = [sv_form_factor_intensity_weight(q, controls) for q in qs]
    out = Matrix{Float64}(I)
    for iq in 1:length(qs)
        out[:, iq] .*= weights[iq]
    end
    return out, weights, amps, qmag
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
    # P1, one-site triangular net in the same 120-degree convention used by
    # the analytical J1/J2 form factors. Because this is P1, do not assume
    # Sunny will generate triangular shell partners from a single representative
    # bond. The explicit shell offsets below are the source of truth.
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    return Crystal(latvecs, [[0.0, 0.0, 0.0]], 1; types=["Yb"])
end

# Positive representatives of the analytical triangular-lattice exchange shells
# for the 120-degree basis used above. Hermitian counterparts are generated by
# Sunny/set_exchange! or set_exchange_at! internally.
#
# These offsets produce the analytical form factors
#   Δ1 = 6 - 2[cos(2πH) + cos(2πK) + cos(2π(H+K))]
#   Δ2 = 6 - 2[cos(2π(H-K)) + cos(2π(2H+K)) + cos(2π(H+2K))]
# so that Δ1(K)=9, Δ2(K)=0 and Δ1(M)=8, Δ2(M)=8.
function sv_j1_shell_offsets()
    return ([1, 0, 0], [0, 1, 0], [1, 1, 0])
end

function sv_j2_shell_offsets()
    return ([1, -1, 0], [2, 1, 0], [1, 2, 0])
end

function sv_offset_string(offsets)
    return join(["(" * join(string.(o), ",") * ")" for o in offsets], "; ")
end

function sv_delta_from_shell_offsets(q::AbstractVector, offsets)
    length(q) == 3 || error("q must have length 3")
    H, K, L = Float64(q[1]), Float64(q[2]), Float64(q[3])
    acc = 0.0
    for o in offsets
        phase = 2π * (H*Float64(o[1]) + K*Float64(o[2]) + L*Float64(o[3]))
        acc += 2.0 * (1.0 - cos(phase))
    end
    return acc
end

function sv_analytical_delta1(q::AbstractVector)
    H, K = Float64(q[1]), Float64(q[2])
    return 6.0 - 2.0 * (cos(2π*H) + cos(2π*K) + cos(2π*(H+K)))
end

function sv_analytical_delta2(q::AbstractVector)
    H, K = Float64(q[1]), Float64(q[2])
    return 6.0 - 2.0 * (cos(2π*(H-K)) + cos(2π*(2H+K)) + cos(2π*(H+2K)))
end

function sv_exchange_geometry_sanity_rows()
    points = [
        ("Γ", [0.0, 0.0, 0.0]),
        ("K_1over3_1over3", [1/3, 1/3, 0.0]),
        ("M_0p5_0_0", [0.5, 0.0, 0.0]),
        ("Γ_equiv_1_0_0", [1.0, 0.0, 0.0]),
    ]
    return [(;
        point = name,
        H = q[1], K = q[2], L = q[3],
        delta1_analytical = sv_analytical_delta1(q),
        delta1_bonds = sv_delta_from_shell_offsets(q, sv_j1_shell_offsets()),
        delta2_analytical = sv_analytical_delta2(q),
        delta2_bonds = sv_delta_from_shell_offsets(q, sv_j2_shell_offsets()),
    ) for (name, q) in points]
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
        # Explicit analytical exchange shells. The Crystal is P1, so do not
        # rely on one representative bond to generate the triangular shell.
        for off in sv_j1_shell_offsets()
            set_exchange!(sys, params.J1_meV, Bond(1, 1, off))
        end
        for off in sv_j2_shell_offsets()
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
        for off in sv_j1_shell_offsets()
            for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
                set_exchange_at!(sys, params.J1_meV * (1 + params.sigma_J * randn(rng)), s1, s2; offset=o)
            end
        end
        for off in sv_j2_shell_offsets()
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
    # Nominal scan coordinates from the CNCS 1D scan file.  These are carried
    # through so Sunny can optionally evaluate finite-Q-averaged spectra around
    # the same measured Q used by the experimental cut rather than only the
    # qtag lookup center.
    H::Vector{Float64}
    K::Vector{Float64}
    L::Vector{Float64}
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
        copy(raw.H),
        copy(raw.K),
        copy(raw.L),
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

# Generic bin-edge helper used for path-coordinate histogram sampling too.
function sv_coordinate_bin_edges(x::AbstractVector)
    n = length(x)
    n >= 2 || error("Need at least two coordinate points to infer bin edges")
    xs = Float64.(x)
    edges = zeros(Float64, n + 1)
    for i in 2:n
        edges[i] = 0.5 * (xs[i-1] + xs[i])
    end
    edges[1] = xs[1] - 0.5 * (xs[2] - xs[1])
    edges[end] = xs[end] + 0.5 * (xs[end] - xs[end-1])
    return edges
end

function sv_lookup_nested_dict(d::Dict, keys::Vector{String}, default)
    cur = d
    for (i, key) in enumerate(keys)
        if !(cur isa Dict) || !haskey(cur, key)
            return default
        end
        if i == length(keys)
            return cur[key]
        end
        cur = cur[key]
    end
    return default
end

function sv_energy_resolution_controls(controls::Dict; section::AbstractString="kpm_2d")
    sec = get(controls, section, Dict{String,Any}())
    er = sec isa Dict ? get(sec, "energy_resolution", Dict{String,Any}()) : Dict{String,Any}()
    er isa Dict || (er = Dict{String,Any}())

    enabled = Bool(get(er, "enabled", false))
    mode = Symbol(get(er, "mode", "none"))
    subtract_kpm_kernel = Bool(get(er, "subtract_kpm_kernel", true))
    kernel_fwhm = Float64(get(get(controls, "kpm", Dict{String,Any}()), "kernel_fwhm_meV", 0.0))
    min_sigma = Float64(get(er, "min_sigma_meV", 1e-6))

    # First-pass CNCS Ei=4.65 meV tabulation used only for post-convolution.
    # The values can be replaced by the analytical-model resolution function later.
    table = if haskey(er, "fwhm_table_meV")
        [Tuple(Float64.(row)) for row in er["fwhm_table_meV"]]
    else
        [(0.5, 0.155), (1.0, 0.136), (1.5, 0.116), (2.0, 0.098), (2.5, 0.083), (3.0, 0.070), (4.0, 0.055)]
    end
    sort!(table; by=x->x[1])
    constant_fwhm = Float64(get(er, "constant_fwhm_meV", kernel_fwhm))
    return (; enabled, mode, subtract_kpm_kernel, kernel_fwhm, min_sigma, table, constant_fwhm)
end

function sv_linear_interp_table(x::Real, table::Vector{Tuple{Float64,Float64}})
    xx = Float64(x)
    isempty(table) && return NaN
    xx <= table[1][1] && return table[1][2]
    xx >= table[end][1] && return table[end][2]
    for i in 1:(length(table)-1)
        x0, y0 = table[i]
        x1, y1 = table[i+1]
        if x0 <= xx <= x1
            t = (xx - x0) / (x1 - x0)
            return (1 - t) * y0 + t * y1
        end
    end
    return table[end][2]
end

function sv_energy_resolution_sigma_meV(E::Real, er)
    if !er.enabled || er.mode in (:none, :off, :disabled)
        return 0.0
    elseif er.mode in (:constant_fwhm, :constant)
        fwhm_target = er.constant_fwhm
    elseif er.mode in (:tabulated_fwhm, :cncs_tabulated, :analytical_like)
        fwhm_target = sv_linear_interp_table(E, er.table)
    else
        error("Unsupported energy-resolution mode $(er.mode). Use none, constant_fwhm, or tabulated_fwhm.")
    end

    sigma_target = fwhm_target / 2.35482004503
    if er.subtract_kpm_kernel
        sigma_kpm = er.kernel_fwhm / 2.35482004503
        sigma_target = sqrt(max(0.0, sigma_target^2 - sigma_kpm^2))
    end
    return max(er.min_sigma, sigma_target)
end

function sv_post_deposit_energy_resolution(E_model::AbstractVector{<:Real}, I_E_by_q::AbstractMatrix, E_target::AbstractVector{<:Real}, er)
    # Deposit the native Sunny/KPM spectrum onto the displayed experimental
    # energy-bin centers using a Gaussian kernel.  This is closer to the
    # analytical event-histogrammer than direct interpolation.  The absolute
    # normalization is not sacred here because the neutron scale is fitted; the
    # important part is the line shape/width.
    if !er.enabled || er.mode in (:none, :off, :disabled)
        return nothing
    end

    Em = Float64.(E_model)
    Et = Float64.(E_target)
    nE, nq = size(I_E_by_q)
    nE == length(Em) || error("Energy axis length mismatch in sv_post_deposit_energy_resolution")
    out = zeros(Float64, length(Et), nq)

    for (im, Etrue) in enumerate(Em)
        sig = sv_energy_resolution_sigma_meV(Etrue, er)
        if !(isfinite(sig) && sig > 0)
            # Nearest-bin fallback.
            j = argmin(abs.(Et .- Etrue))
            out[j, :] .+= I_E_by_q[im, :]
            continue
        end
        w = exp.(-0.5 .* ((Et .- Etrue) ./ sig).^2)
        sw = sum(w)
        if isfinite(sw) && sw > 0
            w ./= sw
            for j in eachindex(Et)
                out[j, :] .+= w[j] .* I_E_by_q[im, :]
            end
        end
    end
    return out
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

function sv_model_to_experimental_energy_grid_resolved(
    E_model::AbstractVector,
    I_model::AbstractVector,
    E_exp::AbstractVector,
    controls::Dict;
    section::AbstractString="kpm",
    mode::Symbol=:bin_average,
)
    er = sv_energy_resolution_controls(controls; section)
    M = reshape(Float64.(I_model), length(I_model), 1)
    deposited = sv_post_deposit_energy_resolution(E_model, M, E_exp, er)
    if deposited !== nothing
        return vec(deposited[:, 1])
    end
    return sv_model_to_experimental_energy_grid(E_model, I_model, E_exp; mode=mode)
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
        println(io, "path,Ei_meV,temperature_K,field_T,qtag,n_E,E_min_meV,E_max_meV,I_min,I_max,rawI_min,rawI_max,bg_min,bg_max,data_mode,background_mean,H_mean,K_mean,L_mean,H_span,K_span,L_span")
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
                mean(c.H),
                mean(c.K),
                mean(c.L),
                maximum(c.H) - minimum(c.H),
                maximum(c.K) - minimum(c.K),
                maximum(c.L) - minimum(c.L),
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
    # qtags are expressed in the same analytical / Mantid-style HKL convention
    # used throughout the co-fit scripts.  Do not apply a K -> -K transform here.
    if qtag in ("0_0_0", "Gamma", "gamma", "Γ")
        return [0.0, 0.0, 0.0]
    elseif qtag == "0_1_0"
        return [0.0, 1.0, 0.0]
    elseif qtag in ("1_0_0", "Gamma_equiv", "gamma_equiv")
        return [1.0, 0.0, 0.0]
    elseif qtag == "0p33_0p33_0"
        return [1/3, 1/3, 0.0]
    elseif qtag == "0p5_0_0"
        return [0.5, 0.0, 0.0]
    else
        error("Unknown qtag $qtag")
    end
end

function sv_qtag_label(qtag::AbstractString)
    if qtag in ("0_0_0", "0_1_0", "1_0_0", "Gamma", "gamma", "Γ", "Gamma_equiv", "gamma_equiv")
        return "Γ"
    elseif qtag == "0p33_0p33_0"
        return "K"
    elseif qtag == "0p5_0_0"
        return "M"
    else
        return String(qtag)
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

# -----------------------------------------------------------------------------
# 1D Sunny KPM finite-Q averaging / experimental-cut histogram approximation
# -----------------------------------------------------------------------------

function sv_kpm_1d_q_averaging_controls(controls::Dict)
    kc = get(controls, "kpm", Dict{String,Any}())
    qa = get(kc, "q_averaging", Dict{String,Any}())
    qa isa Dict || (qa = Dict{String,Any}())
    enabled = Bool(get(qa, "enabled", false))
    mode = Symbol(get(qa, "mode", "gaussian_grid"))
    n_h = Int(get(qa, "n_h", get(qa, "n_H", 3)))
    n_k = Int(get(qa, "n_k", get(qa, "n_K", 3)))
    n_l = Int(get(qa, "n_l", get(qa, "n_L", 1)))
    sigma_H = Float64(get(qa, "sigma_H_rlu", 0.0))
    sigma_K = Float64(get(qa, "sigma_K_rlu", 0.0))
    sigma_L = Float64(get(qa, "sigma_L_rlu", 0.0))
    grid_nsigma = Float64(get(qa, "grid_nsigma", 1.5))
    return (; enabled, mode, n_h, n_k, n_l, sigma_H, sigma_K, sigma_L, grid_nsigma)
end

function sv_kpm_1d_experimental_histogram_controls(controls::Dict)
    kc = get(controls, "kpm", Dict{String,Any}())
    eh = get(kc, "experimental_histogram", Dict{String,Any}())
    eh isa Dict || (eh = Dict{String,Any}())
    enabled = Bool(get(eh, "enabled", false))
    mode = Symbol(get(eh, "mode", "cut_center_plus_resolution_grid"))
    use_scan_q_columns = Bool(get(eh, "use_scan_q_columns", true))

    # Event-style / analytical-cut-volume approximation.  This samples measured
    # Q inside the same nominal 1D cut boxes as the analytical histogrammer, then
    # applies the [kpm.q_averaging] momentum-resolution offsets around each
    # measured Q.  Keep the default grids small: each Sunny KPM call evaluates all
    # q points in one batch, but the cost still scales with n_measured*n_qavg.
    n_measured_h = Int(get(eh, "n_measured_h", get(eh, "n_H", 3)))
    n_measured_k = Int(get(eh, "n_measured_k", get(eh, "n_K", 3)))
    n_measured_l = Int(get(eh, "n_measured_l", get(eh, "n_L", 1)))
    measured_grid_mode = Symbol(get(eh, "measured_grid_mode", "uniform_grid"))
    max_q_samples = Int(get(eh, "max_q_samples", 5000))

    return (; enabled, mode, use_scan_q_columns, n_measured_h, n_measured_k,
        n_measured_l, measured_grid_mode, max_q_samples)
end

function sv_kpm_1d_q_average_offsets(controls::Dict)
    ctl = sv_kpm_1d_q_averaging_controls(controls)
    if !ctl.enabled
        return (; enabled=false, mode=ctl.mode, offsets=[[0.0, 0.0, 0.0]], weights=[1.0],
            n_samples=1, sigma_H=ctl.sigma_H, sigma_K=ctl.sigma_K, sigma_L=ctl.sigma_L,
            n_h=1, n_k=1, n_l=1, grid_nsigma=ctl.grid_nsigma)
    end
    ctl.mode == :gaussian_grid || error("Unsupported [kpm.q_averaging].mode=$(ctl.mode). Currently use gaussian_grid.")
    hs, wh = sv_gaussian_grid_axis(ctl.n_h, ctl.sigma_H, ctl.grid_nsigma)
    ks, wk = sv_gaussian_grid_axis(ctl.n_k, ctl.sigma_K, ctl.grid_nsigma)
    ls, wl = sv_gaussian_grid_axis(ctl.n_l, ctl.sigma_L, ctl.grid_nsigma)
    offsets = Vector{Vector{Float64}}()
    weights = Float64[]
    for (ih, dh) in enumerate(hs), (ik, dk) in enumerate(ks), (il, dl) in enumerate(ls)
        push!(offsets, [Float64(dh), Float64(dk), Float64(dl)])
        push!(weights, wh[ih] * wk[ik] * wl[il])
    end
    sw = sum(weights)
    if !(isfinite(sw) && sw > 0)
        weights .= 1.0 / length(weights)
    else
        weights ./= sw
    end
    return (; enabled=true, mode=ctl.mode, offsets, weights, n_samples=length(weights),
        sigma_H=ctl.sigma_H, sigma_K=ctl.sigma_K, sigma_L=ctl.sigma_L,
        n_h=length(hs), n_k=length(ks), n_l=length(ls), grid_nsigma=ctl.grid_nsigma)
end

function sv_kpm_1d_q_center(cut::SVNeutronCut1D, controls::Dict)
    eh = sv_kpm_1d_experimental_histogram_controls(controls)
    if eh.enabled && eh.use_scan_q_columns && !isempty(cut.H)
        return [mean(cut.H), mean(cut.K), mean(cut.L)]
    else
        return sv_qtag_to_q(cut.qtag)
    end
end

function sv_1d_analytical_cut_ranges(qtag::AbstractString)
    # Same nominal 1D cut boxes used by the analytical histogrammer
    # default_cuts_1d() in the legacy co-fit script.  These ranges describe the
    # measured Q volume being integrated before momentum-resolution convolution.
    if qtag == "0_1_0"
        return (;
            H=(-0.1, 0.1),
            K=(0.9, 1.1),
            L=(-0.3, 0.3),
            label="Gamma_cut_center_0_1_0",
        )
    elseif qtag == "0p33_0p33_0"
        return (;
            H=(0.23, 0.43),
            K=(0.23, 0.43),
            L=(-0.3, 0.3),
            label="M_label_cut_center_1over3_1over3_0",
        )
    elseif qtag == "0p5_0_0"
        return (;
            H=(0.4, 0.6),
            K=(-0.1, 0.1),
            L=(-0.3, 0.3),
            label="K_label_cut_center_1over2_0_0",
        )
    else
        error("No analytical 1D cut-volume ranges defined for qtag=$qtag")
    end
end

function sv_uniform_grid_axis_1d(n::Integer, range_tuple::Tuple{Float64,Float64})
    n = Int(n)
    lo, hi = range_tuple
    if n <= 1 || hi == lo
        return ([0.5 * (lo + hi)], [1.0])
    end
    # Midpoint rule avoids placing samples exactly on sharp cut edges.
    dx = (hi - lo) / n
    xs = [lo + (i - 0.5) * dx for i in 1:n]
    ws = fill(1.0 / n, n)
    return (xs, ws)
end

function sv_kpm_1d_q_sampler(cut::SVNeutronCut1D, controls::Dict)
    eh = sv_kpm_1d_experimental_histogram_controls(controls)
    qavg = sv_kpm_1d_q_average_offsets(controls)

    measured_qs = Vector{Vector{Float64}}()
    measured_weights = Float64[]
    cut_label = "center"
    measured_mode = :center

    if eh.enabled && eh.mode in (:analytical_cut_volume_grid, :event_style_cut_volume, :cut_volume_plus_resolution_grid)
        eh.measured_grid_mode == :uniform_grid || error("Unsupported [kpm.experimental_histogram].measured_grid_mode=$(eh.measured_grid_mode). Use uniform_grid.")
        ranges = sv_1d_analytical_cut_ranges(cut.qtag)
        hs, wh = sv_uniform_grid_axis_1d(eh.n_measured_h, ranges.H)
        ks, wk = sv_uniform_grid_axis_1d(eh.n_measured_k, ranges.K)
        ls, wl = sv_uniform_grid_axis_1d(eh.n_measured_l, ranges.L)
        for (ih, Hm) in enumerate(hs), (ik, Km) in enumerate(ks), (il, Lm) in enumerate(ls)
            push!(measured_qs, [Float64(Hm), Float64(Km), Float64(Lm)])
            push!(measured_weights, wh[ih] * wk[ik] * wl[il])
        end
        cut_label = ranges.label
        measured_mode = :analytical_cut_volume_grid
    else
        push!(measured_qs, sv_kpm_1d_q_center(cut, controls))
        push!(measured_weights, 1.0)
        measured_mode = eh.enabled ? eh.mode : :qtag_center
    end

    qs = Vector{Vector{Float64}}()
    weights = Float64[]
    for (im, qm) in enumerate(measured_qs), (ioff, off) in enumerate(qavg.offsets)
        push!(qs, Float64.(qm) .+ off)
        push!(weights, measured_weights[im] * qavg.weights[ioff])
    end

    sw = sum(weights)
    if !(isfinite(sw) && sw > 0)
        weights .= 1.0 / length(weights)
    else
        weights ./= sw
    end

    if length(qs) > eh.max_q_samples
        error("1D Sunny event-style q sampler requested $(length(qs)) q points, exceeding max_q_samples=$(eh.max_q_samples). Reduce n_measured_* or q_averaging grid, or raise max_q_samples.")
    end

    q_center = [sum(weights[i] * qs[i][j] for i in eachindex(qs)) for j in 1:3]
    return (;
        enabled = (eh.enabled || qavg.enabled),
        mode = eh.enabled ? eh.mode : qavg.mode,
        measured_mode,
        cut_label,
        qs,
        weights,
        q_center,
        n_samples = length(qs),
        n_measured = length(measured_qs),
        n_resolution = qavg.n_samples,
        q_average_enabled = qavg.enabled,
        experimental_histogram_enabled = eh.enabled,
        use_scan_q_columns = eh.use_scan_q_columns,
        sigma_H = qavg.sigma_H,
        sigma_K = qavg.sigma_K,
        sigma_L = qavg.sigma_L,
        n_h = qavg.n_h,
        n_k = qavg.n_k,
        n_l = qavg.n_l,
        grid_nsigma = qavg.grid_nsigma,
        measured_n_h = eh.n_measured_h,
        measured_n_k = eh.n_measured_k,
        measured_n_l = eh.n_measured_l,
        max_q_samples = eh.max_q_samples,
    )
end

function sv_kpm_1d_average_qsampled_intensity(I::AbstractMatrix, sampler)
    nE, nq = size(I)
    nq == sampler.n_samples || error("Expected $(sampler.n_samples) 1D q-sampled columns but got $nq")
    out = zeros(Float64, nE)
    for iq in 1:nq
        out .+= sampler.weights[iq] .* I[:, iq]
    end
    return out
end

function sv_kpm_component_spectrum(params, controls::Dict; component::Symbol, field_T::Real, qtag::AbstractString, cut=nothing)
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

    sampler = if cut === nothing
        qavg0 = sv_kpm_1d_q_average_offsets(controls)
        q0 = sv_qtag_to_q(qtag)
        qs0 = [Float64.(q0) .+ off for off in qavg0.offsets]
        (; enabled=qavg0.enabled, mode=qavg0.mode, measured_mode=:qtag_center,
            cut_label="qtag_center", qs=qs0, weights=qavg0.weights, q_center=Float64.(q0),
            n_samples=qavg0.n_samples, n_measured=1, n_resolution=qavg0.n_samples,
            q_average_enabled=qavg0.enabled, experimental_histogram_enabled=false,
            use_scan_q_columns=false, sigma_H=qavg0.sigma_H, sigma_K=qavg0.sigma_K, sigma_L=qavg0.sigma_L,
            n_h=qavg0.n_h, n_k=qavg0.n_k, n_l=qavg0.n_l, grid_nsigma=qavg0.grid_nsigma,
            measured_n_h=1, measured_n_k=1, measured_n_l=1, max_q_samples=typemax(Int))
    else
        sv_kpm_1d_q_sampler(cut, controls)
    end
    qs = sampler.qs

    energies = collect(range(Float64(kc["energy_min_meV"]), Float64(kc["energy_max_meV"]); length=Int(kc["n_energy"])))
    kernel = gaussian(fwhm=Float64(kc["kernel_fwhm_meV"]))
    swt = SpinWaveTheoryKPM(sys; measure=sv_sunny_measure(sys, controls), tol=Float64(kc["tol"]))
    res = intensities(swt, qs; energies, kernel)
    raw = sv_try_extract_sunny_intensity(res)

    I0 = sv_orient_sunny_intensity_matrix(raw, length(energies), length(qs))
    Ipost, form_factor_weight, form_factor_amplitude, qmag_Ainv = sv_apply_form_factor_to_intensity(I0, qs, controls)
    Iavg = sv_kpm_1d_average_qsampled_intensity(Ipost, sampler)

    return (;
        energy_meV = energies,
        intensity = Iavg,
        intensity_qsampled = Ipost,
        intensity_no_form_factor = I0,
        result = res,
        q = Float64.(sampler.q_center),
        qs = qs,
        q_average = sampler,
        q_average_enabled = sampler.q_average_enabled,
        q_samples = sampler.n_samples,
        q_measured_samples = sampler.n_measured,
        q_resolution_samples = sampler.n_resolution,
        form_factor_weight = sum(sampler.weights .* form_factor_weight),
        form_factor_amplitude = sum(sampler.weights .* form_factor_amplitude),
        qmag_Ainv = sum(sampler.weights .* qmag_Ainv),
    )
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

    ffctl = sv_form_factor_controls(controls)
    fflat = sv_form_factor_lattice_controls(controls)
    er1d = sv_energy_resolution_controls(controls; section="kpm")
    qavg1d = sv_kpm_1d_q_average_offsets(controls)
    eh1d = sv_kpm_1d_experimental_histogram_controls(controls)
    @info "KPM comparison convention" scale_mode scale_scope fallback_neutron_scale flat_to_dispersive_fraction flat_weight histogram_mode=hist_mode energy_resolution_enabled=er1d.enabled energy_resolution_mode=er1d.mode energy_resolution_subtract_kpm_kernel=er1d.subtract_kpm_kernel q_average_enabled=qavg1d.enabled q_average_samples=qavg1d.n_samples sigma_H_rlu=qavg1d.sigma_H sigma_K_rlu=qavg1d.sigma_K sigma_L_rlu=qavg1d.sigma_L experimental_histogram_enabled=eh1d.enabled experimental_histogram_mode=eh1d.mode use_scan_q_columns=eh1d.use_scan_q_columns sunny_transverse_gxy dims=kpm_sizectl.dims repeat_factor=kpm_sizectl.repeat_factor system_size=kpm_sizectl.system_size J1_bonds=sv_offset_string(sv_j1_shell_offsets()) J2_bonds=sv_offset_string(sv_j2_shell_offsets()) data_cuts=length(cuts) form_factor_enabled=ffctl.enabled form_factor_source=ffctl.source form_factor_ion=ffctl.ion form_factor_candidates=ffctl.candidate_ions form_factor_manual_include_j2=ffctl.include_j2 form_factor_manual_apply_as=ffctl.apply_as form_factor_a_A=fflat.a_A form_factor_c_A=fflat.c_A

    # First compute all unscaled Sunny spectra on the experimental energy grids.
    # The neutron intensity scale can then be either fixed, fitted globally, or
    # fitted per cut without rerunning Sunny/KPM.
    cut_results = NamedTuple[]
    for B in fields
        for qtag in qtags
            cut = sv_find_cut(cuts, B, qtag)
            cut === nothing && error("No experimental 1D cut found for B=$B T qtag=$qtag. Check kpm Ei/T/field/qtag controls.")

            @info "Computing Sunny KPM cut and histogramming to experiment" B_T=B qtag=qtag nE_exp=length(cut.energy_meV)
            disp = sv_kpm_component_spectrum(params, controls; component=:dispersive, field_T=B, qtag, cut=cut)
            flat = sv_kpm_component_spectrum(params, controls; component=:flat, field_T=B, qtag, cut=cut)

            Idisp_grid = sv_model_to_experimental_energy_grid_resolved(disp.energy_meV, disp.intensity, cut.energy_meV, controls; section="kpm", mode=hist_mode)
            Iflat_grid = sv_model_to_experimental_energy_grid_resolved(flat.energy_meV, flat.intensity, cut.energy_meV, controls; section="kpm", mode=hist_mode)
            Itotal_unscaled = Idisp_grid .+ flat_weight .* Iflat_grid

            push!(cut_results, (;
                B_T = B,
                qtag = qtag,
                cut = cut,
                Idisp_grid = Idisp_grid,
                Iflat_grid = Iflat_grid,
                Itotal_unscaled = Itotal_unscaled,
                q_center = disp.q,
                qmag_Ainv = disp.qmag_Ainv,
                form_factor_weight = disp.form_factor_weight,
                form_factor_amplitude = disp.form_factor_amplitude,
                q_average_enabled = disp.q_average_enabled,
                q_samples = disp.q_samples,
                q_measured_samples = disp.q_measured_samples,
                q_resolution_samples = disp.q_resolution_samples,
                q_sigma_H = disp.q_average.sigma_H,
                q_sigma_K = disp.q_average.sigma_K,
                q_sigma_L = disp.q_average.sigma_L,
                experimental_histogram_enabled = eh1d.enabled,
                experimental_histogram_mode = eh1d.mode,
                experimental_histogram_measured_mode = disp.q_average.measured_mode,
                experimental_histogram_cut_label = disp.q_average.cut_label,
                use_scan_q_columns = eh1d.use_scan_q_columns,
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
            "energy_meV,I_exp,Ierr_exp,I_total_scaled,I_disp_scaled,I_flat_scaled,I_total_unscaled,I_disp_unscaled,I_flat_unweighted,residual,raw_intensity_exp,background,qx,qy,qz,Q_Ainv_center,form_factor_center,form_factor_weight_center,form_factor_enabled,form_factor_source,form_factor_ion,neutron_scale,neutron_scale_mode,neutron_scale_scope,fallback_neutron_global_scale,flat_weight,flat_to_dispersive_fraction,gperp_ratio,sunny_dims_x,sunny_dims_y,sunny_dims_z,sunny_repeat_x,sunny_repeat_y,sunny_repeat_z,sunny_system_size_x,sunny_system_size_y,sunny_system_size_z,energy_resolution_enabled,energy_resolution_mode,energy_resolution_subtract_kpm_kernel,experimental_histogram_enabled,experimental_histogram_mode,use_scan_q_columns,q_average_enabled,q_samples,q_measured_samples,q_resolution_samples,q_sigma_H_rlu,q_sigma_K_rlu,q_sigma_L_rlu,experimental_histogram_measured_mode,experimental_histogram_cut_label",
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
            fill(r.q_center[1], length(cut.energy_meV)),
            fill(r.q_center[2], length(cut.energy_meV)),
            fill(r.q_center[3], length(cut.energy_meV)),
            fill(r.qmag_Ainv, length(cut.energy_meV)),
            fill(r.form_factor_amplitude, length(cut.energy_meV)),
            fill(r.form_factor_weight, length(cut.energy_meV)),
            fill(ffctl.enabled, length(cut.energy_meV)),
            fill(String(ffctl.source), length(cut.energy_meV)),
            fill(ffctl.ion, length(cut.energy_meV)),
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
            fill(er1d.enabled, length(cut.energy_meV)),
            fill(String(er1d.mode), length(cut.energy_meV)),
            fill(er1d.subtract_kpm_kernel, length(cut.energy_meV)),
            fill(r.experimental_histogram_enabled, length(cut.energy_meV)),
            fill(String(r.experimental_histogram_mode), length(cut.energy_meV)),
            fill(r.use_scan_q_columns, length(cut.energy_meV)),
            fill(r.q_average_enabled, length(cut.energy_meV)),
            fill(r.q_samples, length(cut.energy_meV)),
            fill(r.q_measured_samples, length(cut.energy_meV)),
            fill(r.q_resolution_samples, length(cut.energy_meV)),
            fill(r.q_sigma_H, length(cut.energy_meV)),
            fill(r.q_sigma_K, length(cut.energy_meV)),
            fill(r.q_sigma_L, length(cut.energy_meV)),
            fill(String(r.experimental_histogram_measured_mode), length(cut.energy_meV)),
            fill(r.experimental_histogram_cut_label, length(cut.energy_meV)),
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


# -----------------------------------------------------------------------------
# Sunny parameter-mapping diagnostics
# -----------------------------------------------------------------------------

function sv_check_sunny_parameter_mapping(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    kpm_sizectl = sv_system_size_controls(controls, "kpm")
    large_sizectl = sv_system_size_controls(controls, "largecell")
    r2 = sv_second_kernel_weight(params, controls)
    flat_weight = sv_flat_neutron_weight(params, controls)
    gxy = sv_sunny_transverse_gxy(controls)
    fallback_neutron_scale = sv_neutron_scale(params, controls)

    println("Sunny / analytical parameter mapping")
    println("------------------------------------")
    println("Parameter TOML: ", path)
    println()
    print_canonical_model_parameters(params)

    println("Hamiltonian parameters passed to Sunny")
    println("  dispersive: gzz=$(params.gzz), J1=$(params.J1_meV) meV, J2=$(params.J2_meV) meV")
    println("              sigma_gzz=$(params.sigma_gzz), sigma_J=$(params.sigma_J)")
    println("              J1 shell offsets = ", sv_offset_string(sv_j1_shell_offsets()))
    println("              J2 shell offsets = ", sv_offset_string(sv_j2_shell_offsets()))
    println("  flat:       gzz2=$(params.gzz2), J1=0, J2=0, sigma_gzz2=$(params.sigma_gzz2)")
    println()
    println("Neutron intensity / mixture conventions")
    println("  second_kernel_relative_intensity / flat_to_dispersive_fraction = ", r2)
    println("  gperp_ratio = ", params.gperp_ratio)
    println("  flat_weight used in neutron KPM = r2 * gperp_ratio^2 = ", flat_weight)
    println("  Sunny transverse gxy gauge = ", gxy)
    println("  reference neutron_global_scale = ", fallback_neutron_scale)
    println()
    println("Sunny finite-size controls")
    println("  KPM       dims=$(kpm_sizectl.dims), repeat_factor=$(kpm_sizectl.repeat_factor), system_size=$(kpm_sizectl.system_size)")
    println("  largecell dims=$(large_sizectl.dims), repeat_factor=$(large_sizectl.repeat_factor), system_size=$(large_sizectl.system_size)")

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    mkpath(out_table_dir)
    out_path = joinpath(out_table_dir, "sunny_parameter_mapping.csv")
    names = String[
        "gzz", "J1_meV", "J2_meV", "sigma_gzz", "sigma_J",
        "gzz2", "sigma_gzz2", "gperp_ratio", "chi_vv_muB_per_T",
        "second_kernel_relative_intensity", "flat_weight_r2_gperp_ratio_sq",
        "neutron_global_scale", "magnetization_global_scale",
        "sunny_transverse_gxy",
        "sunny_J1_shell_offsets", "sunny_J2_shell_offsets",
        "geometry_delta1_K_bonds", "geometry_delta2_K_bonds",
        "geometry_delta1_M_bonds", "geometry_delta2_M_bonds",
        "kpm_dims", "kpm_repeat_factor", "kpm_system_size",
        "largecell_dims", "largecell_repeat_factor", "largecell_system_size",
    ]
    values = Any[
        params.gzz, params.J1_meV, params.J2_meV, params.sigma_gzz, params.sigma_J,
        params.gzz2, params.sigma_gzz2, params.gperp_ratio, params.chi_vv_muB_per_T,
        r2, flat_weight, params.neutron_global_scale, params.magnetization_global_scale,
        gxy,
        sv_offset_string(sv_j1_shell_offsets()), sv_offset_string(sv_j2_shell_offsets()),
        sv_delta_from_shell_offsets([1/3, 1/3, 0.0], sv_j1_shell_offsets()), sv_delta_from_shell_offsets([1/3, 1/3, 0.0], sv_j2_shell_offsets()),
        sv_delta_from_shell_offsets([0.5, 0.0, 0.0], sv_j1_shell_offsets()), sv_delta_from_shell_offsets([0.5, 0.0, 0.0], sv_j2_shell_offsets()),
        string(kpm_sizectl.dims), string(kpm_sizectl.repeat_factor), string(kpm_sizectl.system_size),
        string(large_sizectl.dims), string(large_sizectl.repeat_factor), string(large_sizectl.system_size),
    ]
    categories = String[
        "canonical", "canonical", "canonical", "canonical", "canonical",
        "canonical", "canonical", "canonical", "canonical",
        "mixture", "derived_neutron_weight",
        "extrinsic", "extrinsic", "sunny_neutron_gauge",
        "sunny_exchange_geometry", "sunny_exchange_geometry",
        "sunny_exchange_geometry", "sunny_exchange_geometry",
        "sunny_exchange_geometry", "sunny_exchange_geometry",
        "sunny_finite_size", "sunny_finite_size", "sunny_finite_size",
        "sunny_finite_size", "sunny_finite_size", "sunny_finite_size",
    ]
    sv_write_xy_csv(out_path, "name,value,category", names, values, categories)

    geom_rows = sv_exchange_geometry_sanity_rows()
    geom_path = joinpath(out_table_dir, "sunny_exchange_geometry_sanity.csv")
    sv_write_xy_csv(geom_path,
        "point,H,K,L,delta1_analytical,delta1_bonds,delta2_analytical,delta2_bonds",
        [r.point for r in geom_rows],
        [r.H for r in geom_rows],
        [r.K for r in geom_rows],
        [r.L for r in geom_rows],
        [r.delta1_analytical for r in geom_rows],
        [r.delta1_bonds for r in geom_rows],
        [r.delta2_analytical for r in geom_rows],
        [r.delta2_bonds for r in geom_rows],
    )

    println()
    println("Wrote mapping table:")
    println(out_path)
    println("Wrote exchange-geometry sanity table:")
    println(geom_path)
    return (; params, path, r2, flat_weight, sunny_transverse_gxy=gxy, kpm_sizectl, large_sizectl, out_path, geom_path)
end


# -----------------------------------------------------------------------------
# KPM 2D path-map utilities
# -----------------------------------------------------------------------------

function sv_kpm_2d_controls(controls::Dict)
    return get(controls, "kpm_2d", Dict{String,Any}())
end


function sv_kpm_2d_q_averaging_controls(k2)
    qa = get(k2, "q_averaging", Dict{String,Any}())
    qa isa Dict || (qa = Dict{String,Any}())
    enabled = Bool(get(qa, "enabled", false))
    mode = Symbol(get(qa, "mode", "gaussian_grid"))
    n_h = Int(get(qa, "n_h", get(qa, "n_H", 3)))
    n_k = Int(get(qa, "n_k", get(qa, "n_K", 3)))
    n_l = Int(get(qa, "n_l", get(qa, "n_L", 1)))
    sigma_H = Float64(get(qa, "sigma_H_rlu", 0.0))
    sigma_K = Float64(get(qa, "sigma_K_rlu", 0.0))
    sigma_L = Float64(get(qa, "sigma_L_rlu", 0.0))
    grid_nsigma = Float64(get(qa, "grid_nsigma", 1.5))
    return (; enabled, mode, n_h, n_k, n_l, sigma_H, sigma_K, sigma_L, grid_nsigma)
end

function sv_gaussian_grid_axis(n::Integer, sigma::Real, grid_nsigma::Real)
    n = Int(n)
    sig = Float64(sigma)
    if n <= 1 || sig <= 0
        return ([0.0], [1.0])
    end
    xs = collect(range(-Float64(grid_nsigma)*sig, Float64(grid_nsigma)*sig; length=n))
    ws = exp.(-0.5 .* (xs ./ sig).^2)
    sw = sum(ws)
    if !(isfinite(sw) && sw > 0)
        ws .= 1.0 / length(ws)
    else
        ws ./= sw
    end
    return (xs, collect(ws))
end

function sv_kpm_2d_q_average_offsets(k2)
    ctl = sv_kpm_2d_q_averaging_controls(k2)
    if !ctl.enabled
        return (; enabled=false, mode=ctl.mode, offsets=[[0.0, 0.0, 0.0]], weights=[1.0],
            n_samples=1, sigma_H=ctl.sigma_H, sigma_K=ctl.sigma_K, sigma_L=ctl.sigma_L,
            n_h=1, n_k=1, n_l=1, grid_nsigma=ctl.grid_nsigma)
    end
    ctl.mode == :gaussian_grid || error("Unsupported [kpm_2d.q_averaging].mode=$(ctl.mode). Currently use gaussian_grid.")
    hs, wh = sv_gaussian_grid_axis(ctl.n_h, ctl.sigma_H, ctl.grid_nsigma)
    ks, wk = sv_gaussian_grid_axis(ctl.n_k, ctl.sigma_K, ctl.grid_nsigma)
    ls, wl = sv_gaussian_grid_axis(ctl.n_l, ctl.sigma_L, ctl.grid_nsigma)
    offsets = Vector{Vector{Float64}}()
    weights = Float64[]
    for (ih, dh) in enumerate(hs), (ik, dk) in enumerate(ks), (il, dl) in enumerate(ls)
        push!(offsets, [Float64(dh), Float64(dk), Float64(dl)])
        push!(weights, wh[ih] * wk[ik] * wl[il])
    end
    sw = sum(weights)
    if !(isfinite(sw) && sw > 0)
        weights .= 1.0 / length(weights)
    else
        weights ./= sw
    end
    return (; enabled=true, mode=ctl.mode, offsets, weights, n_samples=length(weights),
        sigma_H=ctl.sigma_H, sigma_K=ctl.sigma_K, sigma_L=ctl.sigma_L,
        n_h=length(hs), n_k=length(ks), n_l=length(ls), grid_nsigma=ctl.grid_nsigma)
end

function sv_expand_qs_for_q_averaging(qs::Vector{Vector{Float64}}, qavg)
    if !qavg.enabled
        return qs
    end
    qflat = Vector{Vector{Float64}}()
    sizehint!(qflat, length(qs) * qavg.n_samples)
    for q in qs
        q0 = Float64.(q)
        for off in qavg.offsets
            push!(qflat, q0 .+ off)
        end
    end
    return qflat
end

function sv_average_qsampled_intensity(Iflat::AbstractMatrix, nq::Integer, qavg)
    nE, ncols = size(Iflat)
    if !qavg.enabled
        ncols == nq || error("Expected $nq q columns but got $ncols")
        return Matrix{Float64}(Iflat)
    end
    ns = qavg.n_samples
    ncols == nq * ns || error("Expected $(nq*ns) q-sampled columns but got $ncols")
    out = zeros(Float64, nE, nq)
    for iq in 1:nq
        base = (iq - 1) * ns
        for is in 1:ns
            out[:, iq] .+= qavg.weights[is] .* Iflat[:, base + is]
        end
    end
    return out
end

function sv_path_subsample_axis(lo::Real, hi::Real, n::Integer)
    n = Int(n)
    if n <= 1 || hi <= lo
        return ([0.5 * (Float64(lo) + Float64(hi))], [1.0])
    end
    xs = collect(range(Float64(lo), Float64(hi); length=n+2))[2:end-1]
    ws = fill(1.0 / length(xs), length(xs))
    return (xs, ws)
end

function sv_kpm_2d_experimental_histogram_controls(k2)
    eh = get(k2, "experimental_histogram", Dict{String,Any}())
    eh isa Dict || (eh = Dict{String,Any}())
    enabled = Bool(get(eh, "enabled", false))
    mode = Symbol(get(eh, "mode", "path_bin_plus_resolution_grid"))
    n_path = Int(get(eh, "n_path", 3))
    sample_path_bin = Bool(get(eh, "sample_path_bin", true))
    return (; enabled, mode, n_path, sample_path_bin)
end

function sv_kpm_2d_q_sampler_from_scan(scan, k2; leg::Integer=1)
    # Builds the Q cloud used by the Sunny 2D data/model comparison.
    # Old behavior: one path-center Q plus optional Gaussian H/K offsets.
    # New histogram-like behavior: sample the displayed path-coordinate bin too,
    # then apply the same Gaussian Q-resolution offsets.
    qavg = sv_kpm_2d_q_average_offsets(k2)
    eh = sv_kpm_2d_experimental_histogram_controls(k2)
    xedges = sv_coordinate_bin_edges(scan.x)

    qs_center = sv_kpm_2d_oldpath_qs_from_x(scan.x; leg=leg)
    qs_flat = Vector{Vector{Float64}}()
    ranges = Vector{UnitRange{Int}}()
    weights_by_center = Vector{Vector{Float64}}()

    for ix in eachindex(scan.x)
        start = length(qs_flat) + 1
        local_weights = Float64[]

        if eh.enabled && eh.sample_path_bin
            us, wu = sv_path_subsample_axis(xedges[ix], xedges[ix+1], eh.n_path)
        else
            us, wu = ([Float64(scan.x[ix])], [1.0])
        end

        for (iu, u) in enumerate(us)
            q_path = sv_kpm_2d_oldpath_qs_from_x([u]; leg=leg)[1]
            for is in eachindex(qavg.offsets)
                push!(qs_flat, Float64.(q_path) .+ qavg.offsets[is])
                push!(local_weights, wu[iu] * qavg.weights[is])
            end
        end

        sw = sum(local_weights)
        if !(isfinite(sw) && sw > 0)
            local_weights .= 1.0 / length(local_weights)
        else
            local_weights ./= sw
        end
        stop = length(qs_flat)
        push!(ranges, start:stop)
        push!(weights_by_center, local_weights)
    end

    n_samples = isempty(ranges) ? 0 : length(first(ranges))
    return (; enabled=(qavg.enabled || eh.enabled), mode=eh.enabled ? eh.mode : qavg.mode,
        qs_center, qs_flat, ranges, weights_by_center, n_samples, n_q_centers=length(qs_center),
        q_average_enabled=qavg.enabled, experimental_histogram_enabled=eh.enabled,
        sigma_H=qavg.sigma_H, sigma_K=qavg.sigma_K, sigma_L=qavg.sigma_L,
        n_h=qavg.n_h, n_k=qavg.n_k, n_l=qavg.n_l, grid_nsigma=qavg.grid_nsigma,
        path_n=eh.enabled && eh.sample_path_bin ? eh.n_path : 1)
end

function sv_average_qsampled_intensity_sampler(Iflat::AbstractMatrix, sampler)
    nE, ncols = size(Iflat)
    length(sampler.qs_flat) == ncols || error("Q sampler length $(length(sampler.qs_flat)) does not match intensity columns $ncols")
    out = zeros(Float64, nE, length(sampler.ranges))
    for iq in eachindex(sampler.ranges)
        cols = sampler.ranges[iq]
        ws = sampler.weights_by_center[iq]
        length(cols) == length(ws) || error("Q sampler weight/range mismatch at iq=$iq")
        for (j, col) in enumerate(cols)
            out[:, iq] .+= ws[j] .* Iflat[:, col]
        end
    end
    return out
end

function sv_kpm_component_spectra_for_qs_sampled(params, controls::Dict; component::Symbol, field_T::Real, sampler)
    spec = sv_kpm_component_spectra_for_qs(params, controls; component, field_T, qs=sampler.qs_flat)
    Iavg = sv_average_qsampled_intensity_sampler(spec.intensity, sampler)
    return (; energy_meV=spec.energy_meV, intensity=Iavg, result=spec.result,
        form_factor_weight=spec.form_factor_weight, form_factor_amplitude=spec.form_factor_amplitude, qmag_Ainv=spec.qmag_Ainv,
        q_averaging=sampler, n_q_centers=length(sampler.qs_center), n_q_evaluated=length(sampler.qs_flat))
end

function sv_kpm_component_spectra_for_qs_qaveraged(params, controls::Dict; component::Symbol, field_T::Real, qs::Vector{Vector{Float64}}, qavg)
    qflat = sv_expand_qs_for_q_averaging(qs, qavg)
    spec = sv_kpm_component_spectra_for_qs(params, controls; component, field_T, qs=qflat)
    Iavg = sv_average_qsampled_intensity(spec.intensity, length(qs), qavg)
    return (; energy_meV=spec.energy_meV, intensity=Iavg, result=spec.result,
        form_factor_weight=spec.form_factor_weight, form_factor_amplitude=spec.form_factor_amplitude, qmag_Ainv=spec.qmag_Ainv,
        q_averaging=qavg, n_q_centers=length(qs), n_q_evaluated=length(qflat))
end

function sv_q_path_from_tags(qtags::Vector{String}; n_per_segment::Integer=41)
    length(qtags) >= 2 || error("Need at least two qtags for a 2D path map")
    n_per_segment >= 2 || error("n_per_segment must be >= 2")

    qs = Vector{Vector{Float64}}()
    x = Float64[]
    tick_positions = Float64[]
    tick_labels = String[]

    xpos = 0.0
    for iseg in 1:(length(qtags)-1)
        q0 = Float64.(sv_qtag_to_q(qtags[iseg]))
        q1 = Float64.(sv_qtag_to_q(qtags[iseg+1]))
        seglen = norm(q1 .- q0)
        if iseg == 1
            push!(tick_positions, xpos)
            push!(tick_labels, sv_qtag_label(qtags[iseg]))
        end
        for j in 1:n_per_segment
            if iseg > 1 && j == 1
                continue
            end
            t = (j - 1) / (n_per_segment - 1)
            q = (1 - t) .* q0 .+ t .* q1
            push!(qs, q)
            push!(x, xpos + t * seglen)
        end
        xpos += seglen
        push!(tick_positions, xpos)
        push!(tick_labels, sv_qtag_label(qtags[iseg+1]))
    end

    return (; qs, x, tick_positions, tick_labels)
end

function sv_orient_sunny_intensity_matrix(raw, nE::Integer, nq::Integer)
    if raw isa AbstractVector
        nq == 1 || error("Sunny returned a vector intensity for nq=$nq; expected a matrix")
        return reshape(Float64.(raw), nE, 1)
    elseif raw isa AbstractMatrix
        r, c = size(raw)
        if r == nE && c == nq
            return Matrix{Float64}(raw)
        elseif r == nq && c == nE
            return permutedims(Matrix{Float64}(raw))
        elseif r == nE && c != nq
            @warn "Sunny intensity matrix has energy rows but unexpected q columns; trimming/padding" size_raw=size(raw) nE nq
            out = fill(NaN, nE, nq)
            n = min(c, nq)
            out[:, 1:n] .= Float64.(raw[:, 1:n])
            return out
        elseif c == nE && r != nq
            @warn "Sunny intensity matrix has energy columns but unexpected q rows; trimming/padding" size_raw=size(raw) nE nq
            tmp = permutedims(Matrix{Float64}(raw))
            out = fill(NaN, nE, nq)
            n = min(size(tmp, 2), nq)
            out[:, 1:n] .= tmp[:, 1:n]
            return out
        else
            v = vec(Float64.(raw))
            if length(v) >= nE * nq
                return reshape(v[1:(nE*nq)], nE, nq)
            end
            error("Could not orient Sunny intensity matrix of size $(size(raw)) for nE=$nE nq=$nq")
        end
    else
        error("Unsupported Sunny intensity container type: $(typeof(raw))")
    end
end

function sv_kpm_component_spectra_for_qs(params, controls::Dict; component::Symbol, field_T::Real, qs::Vector{Vector{Float64}})
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

    energies = collect(range(Float64(kc["energy_min_meV"]), Float64(kc["energy_max_meV"]); length=Int(kc["n_energy"])))
    kernel = gaussian(fwhm=Float64(kc["kernel_fwhm_meV"]))
    swt = SpinWaveTheoryKPM(sys; measure=sv_sunny_measure(sys, controls), tol=Float64(kc["tol"]))
    res = intensities(swt, qs; energies, kernel)
    raw = sv_try_extract_sunny_intensity(res)
    I0 = sv_orient_sunny_intensity_matrix(raw, length(energies), length(qs))
    I, form_factor_weight, form_factor_amplitude, qmag_Ainv = sv_apply_form_factor_to_intensity(I0, qs, controls)
    return (; energy_meV=energies, intensity=I, intensity_no_form_factor=I0, result=res,
        form_factor_weight, form_factor_amplitude, qmag_Ainv)
end

function sv_kpm_2d_neutron_scale(params, controls::Dict)
    k2 = sv_kpm_2d_controls(controls)
    mode = Symbol(get(k2, "neutron_scale_mode", "best_fit"))
    if mode == :best_fit
        return sv_neutron_scale(params, controls)
    elseif mode == :manual
        return Float64(get(k2, "manual_neutron_scale", get(controls["kpm"], "manual_neutron_global_scale", 1.0)))
    else
        @warn "kpm_2d neutron_scale_mode=$mode requires experimental data and is not available for a model-only 2D map; using best_fit/reference scale"
        return sv_neutron_scale(params, controls)
    end
end


# -----------------------------------------------------------------------------
# Sunny KPM 2D data-vs-model comparison helpers
# -----------------------------------------------------------------------------

struct SVScan2DCompare
    file::String
    header::String
    xlabel::String
    x::Vector{Float64}
    e::Vector{Float64}
    z::Matrix{Float64}  # size = (length(x), length(e))
end

function sv_header_line_2d(file::AbstractString)
    open(file, "r") do io
        return strip(readline(io))
    end
end

function sv_path_label_from_header_2d(header::AbstractString)
    m = match(r"Error\s+(.+?)\s+DeltaE", header)
    return m === nothing ? "Path coordinate (rlu)" : "Path coordinate " * strip(m.captures[1]) * " (rlu)"
end

function sv_plot_2d_controls(repo_root::AbstractString)
    path = joinpath(repo_root, "configs", "plot_2d_controls.toml")
    return isfile(path) ? load_toml_config(path) : Dict{String,Any}()
end

function sv_neutron_2d_dir(repo_root::AbstractString, controls::Dict)
    # Prefer a Sunny-specific path if present, otherwise reuse the existing
    # plot_2d_controls.toml path so there is only one source for the experimental
    # 2D data location.
    if haskey(controls, "paths") && haskey(controls["paths"], "neutron_2d_dir")
        return sv_repo_path(repo_root, String(controls["paths"]["neutron_2d_dir"]))
    end
    p2 = sv_plot_2d_controls(repo_root)
    if haskey(p2, "data") && haskey(p2["data"], "neutron_2d_subdir")
        return sv_repo_path(repo_root, String(p2["data"]["neutron_2d_subdir"]))
    end
    return joinpath(repo_root, "data", "neutron", "CNCS_2d_scans")
end

function sv_find_scan_file_2d(repo_root::AbstractString, controls::Dict; field_T::Real, leg::Integer=1)
    k2 = sv_kpm_2d_controls(controls)
    data_dir = sv_neutron_2d_dir(repo_root, controls)
    ei_tag = String(get(k2, "ei_tag", "4p65"))
    temp_tag = String(get(k2, "temperature_tag", "0p07K"))
    ft = round(Int, Float64(field_T))
    pattern = Regex("^yzgo_$(ei_tag)meV_$(temp_tag)_$(ft)T_2d_leg$(leg)_SYM\\.dat\$")
    files = sort(filter(f -> occursin(pattern, basename(f)), readdir(data_dir; join=true)))
    if isempty(files)
        return nothing
    elseif length(files) > 1
        @warn "Multiple 2D files matched; using the first" field_T leg files
        return files[1]
    else
        return files[1]
    end
end

function sv_read_scan2d_compare(file::AbstractString; mask_zero::Bool=false)
    header = sv_header_line_2d(file)
    data = readdlm(file, Float64; comments=true, comment_char='#')
    size(data, 2) >= 4 || error("Expected at least 4 numeric columns in $(file), got $(size(data, 2)).")

    intensity = vec(data[:, 1])
    xcol = vec(data[:, 3])
    ecol = vec(data[:, 4])
    xs = sort(collect(unique(xcol)))
    es = sort(collect(unique(ecol)))
    xindex = Dict(v => i for (i, v) in enumerate(xs))
    eindex = Dict(v => i for (i, v) in enumerate(es))
    z = fill(NaN, length(xs), length(es))
    for r in axes(data, 1)
        val = Float64(intensity[r])
        if mask_zero && iszero(val)
            val = NaN
        end
        z[xindex[xcol[r]], eindex[ecol[r]]] = val
    end
    return SVScan2DCompare(String(file), header, sv_path_label_from_header_2d(header), xs, es, z)
end

function sv_load_2d_scans_for_kpm(repo_root::AbstractString, controls::Dict; fields_T::Vector{Float64}, leg::Integer=1)
    scans = Dict{Float64,SVScan2DCompare}()
    for B in fields_T
        file = sv_find_scan_file_2d(repo_root, controls; field_T=B, leg=leg)
        if file === nothing
            @warn "Missing requested 2D scan" field_T=B leg data_dir=sv_neutron_2d_dir(repo_root, controls)
            continue
        end
        scans[B] = sv_read_scan2d_compare(file; mask_zero=get(sv_kpm_2d_controls(controls), "mask_zero_intensity", false))
    end
    return scans
end

function sv_kpm_2d_oldpath_qs_from_x(x::AbstractVector{<:Real}; leg::Integer=1)
    # Matches the old analytical 2D comparison convention for leg 1:
    #   Q(u) = u*[1,-1/2,0] + 1/2*[0,1,0]
    #        = (u, 1/2 - u/2, 0)
    # so u = -1/3, 0, 1/3, 1 correspond to K1, M1, K, Gamma1.
    # Leg 2 is left as a simple placeholder for future extension.
    qs = Vector{Vector{Float64}}()
    for u0 in x
        u = Float64(u0)
        if leg == 1
            push!(qs, [u, 0.5 - 0.5*u, 0.0])
        else
            error("sv_kpm_2d_oldpath_qs_from_x currently implements leg=1 only")
        end
    end
    return qs
end

function sv_finite_values_2d(z)
    vals = Float64[]
    for v in z
        isfinite(v) && push!(vals, Float64(v))
    end
    return vals
end

function sv_quantile_sorted(vals::Vector{Float64}, q::Real)
    isempty(vals) && return NaN
    xs = sort(vals)
    n = length(xs)
    n == 1 && return xs[1]
    t = clamp(Float64(q), 0.0, 1.0) * (n - 1) + 1
    i = floor(Int, t)
    j = ceil(Int, t)
    i == j && return xs[i]
    return xs[i] + (t - i) * (xs[j] - xs[i])
end

function sv_robust_colorrange_2d(arrays; high_quantile::Real=0.995, low_quantile::Real=0.01, force_lo_zero::Bool=true)
    vals = Float64[]
    for z in arrays
        append!(vals, sv_finite_values_2d(z))
    end
    if isempty(vals)
        return (0.0, 1.0)
    end
    lo = force_lo_zero ? 0.0 : sv_quantile_sorted(vals, low_quantile)
    hi = sv_quantile_sorted(vals, high_quantile)
    if !isfinite(lo) || !isfinite(hi) || hi <= lo
        lo = force_lo_zero ? 0.0 : minimum(vals)
        hi = maximum(vals)
    end
    hi <= lo && (hi = lo + 1.0)
    return (lo, hi)
end

function sv_data_colorrange_2d(scans; emin::Real=0.25, emax::Real=Inf, high_quantile::Real=0.95)
    vals = Float64[]
    for s in scans
        for (ie, E) in enumerate(s.e)
            if Float64(emin) <= E <= Float64(emax)
                for ix in eachindex(s.x)
                    v = s.z[ix, ie]
                    isfinite(v) && push!(vals, Float64(v))
                end
            end
        end
    end
    if isempty(vals)
        for s in scans
            append!(vals, sv_finite_values_2d(s.z))
        end
    end
    return sv_robust_colorrange_2d([reshape(vals, length(vals), 1)]; high_quantile=high_quantile)
end

function sv_interp1_linear(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    xx = Float64.(x)
    yy = Float64.(y)
    q = Float64(xq)
    if q < first(xx) || q > last(xx)
        return NaN
    end
    i = searchsortedlast(xx, q)
    if i <= 0
        return yy[1]
    elseif i >= length(xx)
        return yy[end]
    else
        x0, x1 = xx[i], xx[i+1]
        y0, y1 = yy[i], yy[i+1]
        x1 == x0 && return y0
        t = (q - x0) / (x1 - x0)
        return (1 - t) * y0 + t * y1
    end
end

function sv_model_to_scan_energy_grid(model_E::AbstractVector{<:Real}, model_I_E_by_x::AbstractMatrix, scan::SVScan2DCompare; controls=nothing, section::AbstractString="kpm_2d")
    nE, nx = size(model_I_E_by_x)
    nx == length(scan.x) || error("Model q grid length $(nx) does not match scan x length $(length(scan.x))")

    if controls !== nothing
        er = sv_energy_resolution_controls(controls; section)
        deposited = sv_post_deposit_energy_resolution(model_E, model_I_E_by_x, scan.e, er)
        if deposited !== nothing
            # deposited layout is (nE_target, nx); scan z layout is (nx, nE_target)
            return permutedims(deposited)
        end
    end

    out = fill(NaN, length(scan.x), length(scan.e))
    for ix in eachindex(scan.x)
        col = view(model_I_E_by_x, :, ix)
        for (ie, E) in enumerate(scan.e)
            out[ix, ie] = sv_interp1_linear(model_E, col, E)
        end
    end
    return out
end

function sv_least_squares_scale_2d(scans, model_on_scan_grid; emin::Real=0.25, emax::Real=3.2, nonnegative::Bool=true, fallback::Real=1.0)
    num = 0.0
    den = 0.0
    for (B, scan) in scans
        haskey(model_on_scan_grid, B) || continue
        m = model_on_scan_grid[B]
        for ix in eachindex(scan.x), ie in eachindex(scan.e)
            E = scan.e[ie]
            if Float64(emin) <= E <= Float64(emax)
                d = scan.z[ix, ie]
                y = m[ix, ie]
                if isfinite(d) && isfinite(y)
                    num += d * y
                    den += y * y
                end
            end
        end
    end
    if !(isfinite(den) && den > 0)
        @warn "Degenerate 2D least-squares scale; using fallback" fallback den
        return Float64(fallback)
    end
    s = num / den
    nonnegative && (s = max(0.0, s))
    isfinite(s) || return Float64(fallback)
    return Float64(s)
end

function sv_kpm_2d_plot_guides(k2)
    return Float64.(get(k2, "guide_xs", [-1/3, 0.0, 1/3, 2/3, 1.0]))
end

function sv_kpm_2d_plot_ticks(k2)
    xs = Float64.(get(k2, "xtick_positions", [-1/3, 0.0, 1/3, 1.0]))
    labs = String.(get(k2, "xtick_labels", ["K₁", "M₁", "K", "Γ₁"]))
    return (xs, labs)
end

function sv_run_kpm_2d_data_model_comparison(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    print_canonical_model_parameters(params)
    @info "Sunny validation: KPM 2D data-vs-model comparison" sunny_version=sv_try_pkgversion(Sunny) params_path=path

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    out_fig_dir = sv_repo_path(repo_root, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir); mkpath(out_fig_dir)

    k2 = sv_kpm_2d_controls(controls)
    fields = haskey(k2, "fields_T") ? Float64.(k2["fields_T"]) : Float64.(controls["common"]["fields_T"])
    leg = Int(get(k2, "leg", 1))
    scans = sv_load_2d_scans_for_kpm(repo_root, controls; fields_T=fields, leg=leg)
    isempty(scans) && error("No experimental 2D scans loaded for Sunny KPM 2D comparison")

    flat_to_dispersive_fraction = sv_second_kernel_weight(params, controls)
    flat_weight = sv_flat_neutron_weight(params, controls)
    reference_scale = sv_neutron_scale(params, controls)
    scale_mode = Symbol(get(k2, "neutron_scale_mode", "global_least_squares"))
    sunny_transverse_gxy = sv_sunny_transverse_gxy(controls)
    kpm_sizectl = sv_system_size_controls(controls, "kpm")
    qavg = sv_kpm_2d_q_average_offsets(k2)
    ehctl = sv_kpm_2d_experimental_histogram_controls(k2)
    erctl = sv_energy_resolution_controls(controls; section="kpm_2d")

    ffctl = sv_form_factor_controls(controls)
    fflat = sv_form_factor_lattice_controls(controls)
    @info "KPM 2D data/model convention" fields leg flat_to_dispersive_fraction flat_weight reference_scale scale_mode sunny_transverse_gxy dims=kpm_sizectl.dims repeat_factor=kpm_sizectl.repeat_factor system_size=kpm_sizectl.system_size J1_bonds=sv_offset_string(sv_j1_shell_offsets()) J2_bonds=sv_offset_string(sv_j2_shell_offsets()) q_average_enabled=qavg.enabled q_average_samples=qavg.n_samples sigma_H_rlu=qavg.sigma_H sigma_K_rlu=qavg.sigma_K sigma_L_rlu=qavg.sigma_L experimental_histogram_enabled=ehctl.enabled path_bin_samples=ehctl.n_path energy_resolution_enabled=erctl.enabled energy_resolution_mode=erctl.mode energy_resolution_subtract_kpm_kernel=erctl.subtract_kpm_kernel form_factor_enabled=ffctl.enabled form_factor_source=ffctl.source form_factor_ion=ffctl.ion form_factor_candidates=ffctl.candidate_ions form_factor_manual_include_j2=ffctl.include_j2 form_factor_manual_apply_as=ffctl.apply_as form_factor_a_A=fflat.a_A form_factor_c_A=fflat.c_A

    raw_models = Dict{Float64,Any}()
    model_on_scan_grid = Dict{Float64,Matrix{Float64}}()
    disp_on_scan_grid = Dict{Float64,Matrix{Float64}}()
    flat_on_scan_grid = Dict{Float64,Matrix{Float64}}()
    for B in fields
        haskey(scans, B) || continue
        scan = scans[B]
        sampler = sv_kpm_2d_q_sampler_from_scan(scan, k2; leg=leg)
        qs = sampler.qs_center
        @info "Computing Sunny KPM 2D map on experimental histogram grid" field_T=B leg n_q=length(qs) q_average_enabled=sampler.q_average_enabled experimental_histogram_enabled=sampler.experimental_histogram_enabled q_samples_per_pixel=sampler.n_samples n_q_evaluated=length(sampler.qs_flat) energy_resolution_enabled=erctl.enabled energy_resolution_mode=erctl.mode
        disp = sv_kpm_component_spectra_for_qs_sampled(params, controls; component=:dispersive, field_T=B, sampler=sampler)
        flat = sv_kpm_component_spectra_for_qs_sampled(params, controls; component=:flat, field_T=B, sampler=sampler)
        Itotal_unscaled = disp.intensity .+ flat_weight .* flat.intensity
        raw_models[B] = (; disp, flat, Itotal_unscaled, qs, q_averaging=sampler)
        model_on_scan_grid[B] = sv_model_to_scan_energy_grid(disp.energy_meV, Itotal_unscaled, scan; controls=controls, section="kpm_2d")
        disp_on_scan_grid[B] = sv_model_to_scan_energy_grid(disp.energy_meV, disp.intensity, scan; controls=controls, section="kpm_2d")
        flat_on_scan_grid[B] = sv_model_to_scan_energy_grid(flat.energy_meV, flat.intensity, scan; controls=controls, section="kpm_2d")
    end

    scale_by_field = Dict{Float64,Float64}()
    if scale_mode == :best_fit
        for B in keys(raw_models); scale_by_field[B] = reference_scale; end
    elseif scale_mode == :manual
        s = Float64(get(k2, "manual_neutron_scale", reference_scale))
        for B in keys(raw_models); scale_by_field[B] = s; end
    elseif scale_mode == :global_least_squares
        s = sv_least_squares_scale_2d(scans, model_on_scan_grid;
            emin=Float64(get(k2, "scale_energy_min_meV", 0.25)),
            emax=Float64(get(k2, "scale_energy_max_meV", 3.2)),
            fallback=reference_scale,
        )
        for B in keys(raw_models); scale_by_field[B] = s; end
    elseif scale_mode == :panel_least_squares
        for B in keys(raw_models)
            scale_by_field[B] = sv_least_squares_scale_2d(Dict(B => scans[B]), Dict(B => model_on_scan_grid[B]);
                emin=Float64(get(k2, "scale_energy_min_meV", 0.25)),
                emax=Float64(get(k2, "scale_energy_max_meV", 3.2)),
                fallback=reference_scale,
            )
        end
    else
        error("Unknown [kpm_2d].neutron_scale_mode=$scale_mode. Use best_fit, manual, global_least_squares, or panel_least_squares.")
    end

    # Save one long CSV per field for reproducibility. Model intensities are on
    # the native Sunny KPM energy grid, while the scaling was computed after
    # interpolation to the experimental energy grid.
    csv_paths = Dict{Float64,String}()
    for B in sort(collect(keys(raw_models)))
        m = raw_models[B]
        scan = scans[B]
        scale = scale_by_field[B]
        nE = length(m.disp.energy_meV)
        nq = length(scan.x)
        rows = nE * nq
        q_index = Int[]; path_coordinate = Float64[]; qx=Float64[]; qy=Float64[]; qz=Float64[]; energy_col=Float64[]
        total_scaled=Float64[]; disp_scaled=Float64[]; flat_scaled=Float64[]; total_unscaled=Float64[]; disp_unscaled=Float64[]; flat_unweighted=Float64[]
        for iq in 1:nq
            q = m.qs[iq]
            for iE in 1:nE
                push!(q_index, iq); push!(path_coordinate, scan.x[iq]); push!(qx, q[1]); push!(qy, q[2]); push!(qz, q[3]); push!(energy_col, m.disp.energy_meV[iE])
                push!(total_unscaled, m.Itotal_unscaled[iE, iq])
                push!(disp_unscaled, m.disp.intensity[iE, iq])
                push!(flat_unweighted, m.flat.intensity[iE, iq])
                push!(total_scaled, scale * m.Itotal_unscaled[iE, iq])
                push!(disp_scaled, scale * m.disp.intensity[iE, iq])
                push!(flat_scaled, scale * flat_weight * m.flat.intensity[iE, iq])
            end
        end
        tag = @sprintf("%gT_leg%d", B, leg)
        csv_path = joinpath(out_table_dir, "sunny_kpm_2d_data_model_$(tag).csv")
        sv_write_xy_csv(csv_path,
            "q_index,path_coordinate,qx,qy,qz,Q_Ainv_center,form_factor_center,form_factor_weight_center,energy_meV,I_total_scaled,I_disp_scaled,I_flat_scaled,I_total_unscaled,I_disp_unscaled,I_flat_unweighted,neutron_scale,neutron_scale_mode,flat_weight,flat_to_dispersive_fraction,gperp_ratio,sunny_transverse_gxy,q_average_enabled,q_average_n_samples,q_average_sigma_H_rlu,q_average_sigma_K_rlu,q_average_sigma_L_rlu,form_factor_enabled",
            q_index, path_coordinate, qx, qy, qz,
            vcat([fill(sv_qmag_Ainv(m.qs[iq], controls), nE) for iq in 1:nq]...),
            vcat([fill(sv_form_factor_amplitude(m.qs[iq], controls), nE) for iq in 1:nq]...),
            vcat([fill(sv_form_factor_intensity_weight(m.qs[iq], controls), nE) for iq in 1:nq]...),
            energy_col,
            total_scaled, disp_scaled, flat_scaled, total_unscaled, disp_unscaled, flat_unweighted,
            fill(scale, rows), fill(String(scale_mode), rows), fill(flat_weight, rows), fill(flat_to_dispersive_fraction, rows), fill(params.gperp_ratio, rows), fill(sunny_transverse_gxy, rows),
            fill(qavg.enabled, rows), fill(qavg.n_samples, rows), fill(qavg.sigma_H, rows), fill(qavg.sigma_K, rows), fill(qavg.sigma_L, rows), fill(ffctl.enabled, rows),
        )
        csv_paths[B] = csv_path

        # Also save the model after deposition/interpolation onto the experimental
        # scan grid.  This is the grid used for scaling and plotting when
        # [kpm_2d.energy_resolution] is enabled.
        scan_csv_path = joinpath(out_table_dir, "sunny_kpm_2d_data_model_$(tag)_scan_grid.csv")
        q_index2 = Int[]; x2 = Float64[]; e2 = Float64[]; qx2 = Float64[]; qy2 = Float64[]; qz2 = Float64[]
        itot2 = Float64[]; idisp2 = Float64[]; iflat2 = Float64[]; data2 = Float64[]
        for iq in 1:nq
            q = m.qs[iq]
            for ie in eachindex(scan.e)
                push!(q_index2, iq); push!(x2, scan.x[iq]); push!(e2, scan.e[ie]); push!(qx2, q[1]); push!(qy2, q[2]); push!(qz2, q[3])
                push!(itot2, scale * model_on_scan_grid[B][iq, ie])
                push!(idisp2, scale * disp_on_scan_grid[B][iq, ie])
                push!(iflat2, scale * flat_weight * flat_on_scan_grid[B][iq, ie])
                push!(data2, scan.z[iq, ie])
            end
        end
        nrows2 = length(e2)
        sv_write_xy_csv(scan_csv_path,
            "q_index,path_coordinate,qx,qy,qz,energy_meV,I_exp,I_total_scaled,I_disp_scaled,I_flat_scaled,neutron_scale,neutron_scale_mode,experimental_histogram_enabled,path_bin_samples,q_samples_per_pixel,energy_resolution_enabled,energy_resolution_mode,energy_resolution_subtract_kpm_kernel",
            q_index2, x2, qx2, qy2, qz2, e2, data2, itot2, idisp2, iflat2,
            fill(scale, nrows2), fill(String(scale_mode), nrows2), fill(m.q_averaging.experimental_histogram_enabled, nrows2), fill(m.q_averaging.path_n, nrows2), fill(m.q_averaging.n_samples, nrows2),
            fill(erctl.enabled, nrows2), fill(String(erctl.mode), nrows2), fill(erctl.subtract_kpm_kernel, nrows2))
    end

    # Linear-intensity, 2-row style matching the analytical 2D data/model plot.
    fields_have = sort(collect(keys(raw_models)))
    ncols = length(fields_have)
    fig = Figure(size=(1180, 850), fontsize=16)
    qavg_label = qavg.enabled ? string(", Q-avg ", qavg.n_samples, " offsets") : ", path centers"
    eh_label = ehctl.enabled ? string(", path-bin ×", ehctl.n_path) : ""
    er_label = erctl.enabled ? string(", E-res ", erctl.mode) : ""
    ff_label = ffctl.enabled ? ", Yb³⁺ |f(Q)|²" : ""
    Label(fig[0, 1:(ncols+1)], @sprintf("YZGO 2D data vs Sunny KPM model, leg %d, Ei=4.65 meV, T=0.07 K%s%s%s%s", leg, qavg_label, eh_label, er_label, ff_label), fontsize=21, font=:bold, tellwidth=false)

    energy_ylim = Float64.(get(k2, "energy_ylim_meV", [0.20, 3.20]))
    xlim = haskey(k2, "xlim") ? Float64.(k2["xlim"]) : nothing
    data_cr = sv_data_colorrange_2d([scans[B] for B in fields_have];
        emin=Float64(get(k2, "data_color_energy_min_meV", 0.25)),
        emax=Float64(get(k2, "data_color_energy_max_meV", Inf)),
        high_quantile=Float64(get(k2, "data_clip_high_quantile", 0.95)),
    )
    model_arrays = [scale_by_field[B] .* model_on_scan_grid[B] for B in fields_have]
    model_cr = sv_robust_colorrange_2d(model_arrays; high_quantile=Float64(get(k2, "model_clip_high_quantile", 0.995)))
    cmap = Symbol(get(k2, "colormap", "viridis"))
    guides = sv_kpm_2d_plot_guides(k2)
    xticks = sv_kpm_2d_plot_ticks(k2)

    data_hm = nothing
    model_hm = nothing
    for (icol, B) in enumerate(fields_have)
        scan = scans[B]
        m = raw_models[B]
        zmodel = scale_by_field[B] .* model_on_scan_grid[B]

        axd = Axis(fig[1, icol], title=@sprintf("%g T", B), ylabel=icol == 1 ? "Data\nΔE (meV)" : "", xlabel="")
        data_hm = heatmap!(axd, scan.x, scan.e, scan.z; colormap=cmap, colorrange=data_cr, nan_color=:lightgray)
        xlim === nothing ? xlims!(axd, minimum(scan.x), maximum(scan.x)) : xlims!(axd, xlim[1], xlim[2])
        length(energy_ylim) == 2 && ylims!(axd, energy_ylim[1], energy_ylim[2])
        vlines!(axd, guides; color=(:white, 0.45), linewidth=1)
        axd.xticks = xticks

        axm = Axis(fig[2, icol], title=@sprintf("scale = %.4g", scale_by_field[B]), ylabel=icol == 1 ? "Sunny KPM\nΔE (meV)" : "", xlabel=scan.xlabel)
        model_hm = heatmap!(axm, scan.x, scan.e, zmodel; colormap=cmap, colorrange=model_cr, nan_color=:lightgray)
        xlim === nothing ? xlims!(axm, minimum(scan.x), maximum(scan.x)) : xlims!(axm, xlim[1], xlim[2])
        length(energy_ylim) == 2 && ylims!(axm, energy_ylim[1], energy_ylim[2])
        vlines!(axm, guides; color=(:white, 0.45), linewidth=1)
        axm.xticks = xticks
    end
    data_hm !== nothing && Colorbar(fig[1, ncols+1], data_hm; label="Data intensity")
    model_hm !== nothing && Colorbar(fig[2, ncols+1], model_hm; label="Scaled Sunny KPM intensity")
    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)

    fig_path = joinpath(out_fig_dir, @sprintf("sunny_kpm_2d_data_vs_model_leg%d.png", leg))
    save(fig_path, fig)
    if get(k2, "save_pdf", false)
        try
            save(joinpath(out_fig_dir, @sprintf("sunny_kpm_2d_data_vs_model_leg%d.pdf", leg)), fig)
        catch err
            @warn "Could not save Sunny KPM 2D PDF" exception=(err, catch_backtrace())
        end
    end

    @info "Saved Sunny KPM 2D data/model comparison" fig_path scale_mode scale_by_field csv_paths
    println()
    println("Sunny KPM 2D data-vs-model comparison completed")
    println("  figure: ", fig_path)
    for B in fields_have
        println(@sprintf("  %g T scale = %.6g; CSV = %s", B, scale_by_field[B], csv_paths[B]))
    end
    return (; fig_path, csv_paths, scans, raw_models, scale_by_field, fields_T=fields_have)
end

function sv_run_kpm_2d_path_map(repo_root::AbstractString; controls=sv_load_controls(repo_root))
    (; params, path) = sv_load_params(repo_root, controls)
    print_canonical_model_parameters(params)
    @info "Sunny validation: KPM 2D path map" sunny_version=sv_try_pkgversion(Sunny) params_path=path

    out_table_dir = sv_repo_path(repo_root, controls["paths"]["table_subdir"])
    out_fig_dir = sv_repo_path(repo_root, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir); mkpath(out_fig_dir)

    k2 = sv_kpm_2d_controls(controls)
    B = Float64(get(k2, "field_T", first(Float64.(controls["common"]["fields_T"]))))
    default_2d_qtags = ["0_0_0", "0p33_0p33_0", "0p5_0_0", "1_0_0"]
    qtags = String.(get(k2, "qtags", default_2d_qtags))
    n_per_segment = Int(get(k2, "n_per_segment", 41))
    pathinfo = sv_q_path_from_tags(qtags; n_per_segment=n_per_segment)

    flat_to_dispersive_fraction = sv_second_kernel_weight(params, controls)
    flat_weight = sv_flat_neutron_weight(params, controls)
    neutron_scale = sv_kpm_2d_neutron_scale(params, controls)
    sunny_transverse_gxy = sv_sunny_transverse_gxy(controls)
    kpm_sizectl = sv_system_size_controls(controls, "kpm")
    qavg = sv_kpm_2d_q_average_offsets(k2)

    ffctl = sv_form_factor_controls(controls)
    fflat = sv_form_factor_lattice_controls(controls)
    @info "KPM 2D convention" B_T=B qtags flat_to_dispersive_fraction flat_weight neutron_scale sunny_transverse_gxy dims=kpm_sizectl.dims repeat_factor=kpm_sizectl.repeat_factor system_size=kpm_sizectl.system_size n_q=length(pathinfo.qs) q_average_enabled=qavg.enabled q_average_samples=qavg.n_samples sigma_H_rlu=qavg.sigma_H sigma_K_rlu=qavg.sigma_K sigma_L_rlu=qavg.sigma_L form_factor_enabled=ffctl.enabled form_factor_source=ffctl.source form_factor_ion=ffctl.ion form_factor_candidates=ffctl.candidate_ions form_factor_manual_include_j2=ffctl.include_j2 form_factor_manual_apply_as=ffctl.apply_as form_factor_a_A=fflat.a_A form_factor_c_A=fflat.c_A

    @info "Computing 2D dispersive KPM map" B_T=B n_q=length(pathinfo.qs) n_q_evaluated=length(pathinfo.qs)*qavg.n_samples
    disp = sv_kpm_component_spectra_for_qs_qaveraged(params, controls; component=:dispersive, field_T=B, qs=pathinfo.qs, qavg=qavg)
    @info "Computing 2D flat KPM map" B_T=B n_q=length(pathinfo.qs) n_q_evaluated=length(pathinfo.qs)*qavg.n_samples
    flat = sv_kpm_component_spectra_for_qs_qaveraged(params, controls; component=:flat, field_T=B, qs=pathinfo.qs, qavg=qavg)

    Itotal_unscaled = disp.intensity .+ flat_weight .* flat.intensity
    Itotal_scaled = neutron_scale .* Itotal_unscaled
    Idisp_scaled = neutron_scale .* disp.intensity
    Iflat_scaled = neutron_scale .* flat_weight .* flat.intensity

    basename_tag = get(k2, "output_tag", @sprintf("%gT", B))
    csv_path = joinpath(out_table_dir, "sunny_kpm_2d_path_$(basename_tag).csv")
    nE = length(disp.energy_meV)
    nq = length(pathinfo.qs)
    rows = nE * nq
    q_index = Int[]
    path_coordinate = Float64[]
    qx = Float64[]
    qy = Float64[]
    qz = Float64[]
    energy_col = Float64[]
    total_scaled_col = Float64[]
    disp_scaled_col = Float64[]
    flat_scaled_col = Float64[]
    total_unscaled_col = Float64[]
    disp_unscaled_col = Float64[]
    flat_unweighted_col = Float64[]

    for iq in 1:nq
        q = pathinfo.qs[iq]
        for iE in 1:nE
            push!(q_index, iq)
            push!(path_coordinate, pathinfo.x[iq])
            push!(qx, q[1]); push!(qy, q[2]); push!(qz, q[3])
            push!(energy_col, disp.energy_meV[iE])
            push!(total_scaled_col, Itotal_scaled[iE, iq])
            push!(disp_scaled_col, Idisp_scaled[iE, iq])
            push!(flat_scaled_col, Iflat_scaled[iE, iq])
            push!(total_unscaled_col, Itotal_unscaled[iE, iq])
            push!(disp_unscaled_col, disp.intensity[iE, iq])
            push!(flat_unweighted_col, flat.intensity[iE, iq])
        end
    end

    sv_write_xy_csv(csv_path,
        "q_index,path_coordinate,qx,qy,qz,Q_Ainv_center,form_factor_center,form_factor_weight_center,energy_meV,I_total_scaled,I_disp_scaled,I_flat_scaled,I_total_unscaled,I_disp_unscaled,I_flat_unweighted,neutron_scale,flat_weight,flat_to_dispersive_fraction,gperp_ratio,sunny_transverse_gxy,form_factor_enabled",
        q_index, path_coordinate, qx, qy, qz,
        vcat([fill(sv_qmag_Ainv(pathinfo.qs[iq], controls), nE) for iq in 1:nq]...),
        vcat([fill(sv_form_factor_amplitude(pathinfo.qs[iq], controls), nE) for iq in 1:nq]...),
        vcat([fill(sv_form_factor_intensity_weight(pathinfo.qs[iq], controls), nE) for iq in 1:nq]...),
        energy_col,
        total_scaled_col, disp_scaled_col, flat_scaled_col, total_unscaled_col, disp_unscaled_col, flat_unweighted_col,
        fill(neutron_scale, rows), fill(flat_weight, rows), fill(flat_to_dispersive_fraction, rows), fill(params.gperp_ratio, rows), fill(sunny_transverse_gxy, rows), fill(ffctl.enabled, rows)
    )

    z_mode = Symbol(get(k2, "z_mode", "linear"))
    z_floor = Float64(get(k2, "log10_floor", 1e-12))
    Z = if z_mode == :log10
        log10.(max.(Itotal_scaled, z_floor))
    elseif z_mode == :linear
        Itotal_scaled
    else
        error("Unknown [kpm_2d] z_mode=$z_mode. Use linear or log10.")
    end

    fig = Figure(size=(1100, 760))
    ax = Axis(fig[1, 1],
        xlabel = "Q path",
        ylabel = "Energy transfer (meV)",
        title = @sprintf("Sunny KPM 2D path map, B = %.3g T", B),
        xticks = (pathinfo.tick_positions, pathinfo.tick_labels),
    )
    hm = heatmap!(ax, pathinfo.x, disp.energy_meV, permutedims(Z))
    for xp in pathinfo.tick_positions
        vlines!(ax, [xp], linewidth=1)
    end
    if haskey(k2, "energy_ylim_meV")
        yl = Float64.(k2["energy_ylim_meV"])
        length(yl) == 2 && ylims!(ax, yl[1], yl[2])
    end
    Colorbar(fig[1, 2], hm, label = z_mode == :log10 ? "log10 intensity" : "Intensity (scaled arb.)")
    fig_path = joinpath(out_fig_dir, "sunny_kpm_2d_path_$(basename_tag).png")
    save(fig_path, fig)

    @info "Saved Sunny KPM 2D path map" fig_path csv_path
    println()
    println("Sunny KPM 2D path map completed")
    println("  field_T: ", B)
    println("  qtags:   ", join(qtags, " -> "))
    println("  figure:  ", fig_path)
    println("  CSV:     ", csv_path)
    return (; fig_path, csv_path, field_T=B, qtags, energy_meV=disp.energy_meV, x=pathinfo.x, Itotal_scaled, Idisp_scaled, Iflat_scaled)
end


end # module
