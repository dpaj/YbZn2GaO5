# =============================================================================
# Feature extraction utilities for YbZn2GaO5
# =============================================================================
#
# This file is part of the refactored repo code.  It intentionally extracts
# empirical/interpretable features from the experimental data, rather than
# performing the full microscopic Hamiltonian co-fit.
#
# Current scope:
#   - neutron 1D CNCS scans at multiple fields, fitted independently by field;
#   - magnetization curve, fitted to a compact empirical two-component model;
#   - one combined feature table with block = neutron / magnetization.
#
# Current science-facing conventions:
#   - neutron flat/nondispersive features are allowed to shift and broaden with
#     field;
#   - no legacy sigma_J1/sigma_J2 or gperp/gperp2 terminology is introduced;
#   - magnetization features are empirical descriptors, not final microscopic
#     fit parameters.
#
# =============================================================================

using DelimitedFiles
using Printf
using Statistics
using LinearAlgebra
using Random
using Optim
using CairoMakie

const FE_MU_B_MEV_PER_T = 5.7883818060e-2
const FE_KB_MEV_PER_K   = 8.617333262e-2
const FE_GAUSS_FWHM = 2.3548200450309493

const FE_Q_SHORT = Dict(
    "0_1_0" => "gamma",
    "0p33_0p33_0" => "k",
    "0p5_0_0" => "m",
)

const FE_Q_LABEL = Dict(
    "0_1_0" => "Q = (0, 1, 0)",
    "0p33_0p33_0" => "Q = (1/3, 1/3, 0)",
    "0p5_0_0" => "Q = (1/2, 0, 0)",
)


# Representative in-plane q centers used for co-fit-informed feature initial guesses.
# These match the simple triangular-lattice convention used in the analytical
# polarized-phase model.  They are used only to seed empirical feature fits, not
# to constrain the final extracted features.
const FE_Q_CENTER_HK = Dict(
    "0_1_0" => (0.0, 1.0),
    "0p33_0p33_0" => (1.0/3.0, 1.0/3.0),
    "0p5_0_0" => (0.5, 0.0),
)

function _fe_delta1_triangular(H::Real, K::Real)
    return 6.0 - 2.0 * (
        cos(2.0 * pi * H) +
        cos(2.0 * pi * K) +
        cos(2.0 * pi * (H + K))
    )
end

function _fe_delta2_triangular(H::Real, K::Real)
    return 6.0 - 2.0 * (
        cos(2.0 * pi * (H - K)) +
        cos(2.0 * pi * (2.0 * H + K)) +
        cos(2.0 * pi * (H + 2.0 * K))
    )
end

function _fe_dispersion_guess_from_cofit(qtag::String, B_T::Real, params)
    haskey(FE_Q_CENTER_HK, qtag) || return NaN
    H, K = FE_Q_CENTER_HK[qtag]
    return params.gzz * FE_MU_B_MEV_PER_T * Float64(B_T) -
           0.5 * (params.J1_meV * _fe_delta1_triangular(H, K) +
                  params.J2_meV * _fe_delta2_triangular(H, K))
end

function _fe_dispersion_sigma_guess_from_cofit(qtag::String, B_T::Real, params;
                                               resolution_sigma_meV::Real=0.08)
    haskey(FE_Q_CENTER_HK, qtag) || return 0.20
    H, K = FE_Q_CENTER_HK[qtag]
    d1 = _fe_delta1_triangular(H, K)
    d2 = _fe_delta2_triangular(H, K)
    zeeman_sigma = abs(params.sigma_gzz * FE_MU_B_MEV_PER_T * Float64(B_T))
    # sigma_J is the shared fractional exchange-disorder width in the current model.
    exchange_sigma = abs(params.sigma_J) * sqrt((0.5 * params.J1_meV * d1)^2 +
                                               (0.5 * params.J2_meV * d2)^2)
    return sqrt(zeeman_sigma^2 + exchange_sigma^2 + Float64(resolution_sigma_meV)^2)
end

function _fe_flat_center_guess_from_cofit(B_T::Real, params)
    return params.gzz2 * FE_MU_B_MEV_PER_T * Float64(B_T)
end

function _fe_flat_sigma_guess_from_cofit(B_T::Real, params; resolution_sigma_meV::Real=0.08)
    g_sigma = abs(params.sigma_gzz2 * FE_MU_B_MEV_PER_T * Float64(B_T))
    return sqrt(g_sigma^2 + Float64(resolution_sigma_meV)^2)
end

function _fe_nearest_intensity_at(s, energy_meV::Real)
    isempty(s.energy) && return NaN
    idx = argmin(abs.(s.energy .- Float64(energy_meV)))
    return s.intensity[idx]
end

# -----------------------------------------------------------------------------
# Small optimization helpers
# -----------------------------------------------------------------------------

Base.@kwdef struct FEParamSpec
    name::Symbol
    lo::Float64
    hi::Float64
    initial::Float64
end

function _fe_safe_initial(x::Real, lo::Real, hi::Real)
    width = Float64(hi - lo)
    width > 0 || error("Bad parameter bounds: lo=$lo hi=$hi")
    epsx = max(1e-9 * width, 1e-12)
    return clamp(Float64(x), Float64(lo) + epsx, Float64(hi) - epsx)
end

function _fe_logit(x::Real)
    y = clamp(Float64(x), 1e-12, 1.0 - 1e-12)
    return log(y / (1.0 - y))
end

function _fe_invlogit_stable(u::Real)
    x = Float64(u)
    if x >= 0
        z = exp(-x)
        return 1.0 / (1.0 + z)
    else
        z = exp(x)
        return z / (1.0 + z)
    end
end

function _fe_unconstrained_initial(specs::Vector{FEParamSpec})
    u0 = Float64[]
    for s in specs
        x0 = _fe_safe_initial(s.initial, s.lo, s.hi)
        frac = (x0 - s.lo) / (s.hi - s.lo)
        push!(u0, _fe_logit(frac))
    end
    return u0
end

function _fe_unpack(u::AbstractVector{<:Real}, specs::Vector{FEParamSpec})
    length(u) == length(specs) || error("Parameter vector length mismatch")
    p = Dict{Symbol, Float64}()
    for (ui, s) in zip(u, specs)
        p[s.name] = s.lo + (s.hi - s.lo) * _fe_invlogit_stable(ui)
    end
    return p
end

# -----------------------------------------------------------------------------
# Neutron data loading
# -----------------------------------------------------------------------------

struct FEScan1D
    qtag::String
    field_T::Float64
    temperature_K::Float64
    Ei_meV::Float64
    filename::String
    intensity::Vector{Float64}
    error::Vector{Float64}
    energy::Vector{Float64}
    K::Vector{Float64}
    L::Vector{Float64}
    H::Vector{Float64}
end

_fe_parsefloat_token(s::AbstractString) = parse(Float64, replace(s, "p" => "."))

function fe_parse_neutron_filename(path::AbstractString)
    fname = basename(path)
    m = match(r"^yzgo_(\d+(?:p\d+)?)meV_(\d+(?:p\d+)?)K_(\d+(?:p\d+)?)T_Escan_(.+)_SYM\.dat$", fname)
    m === nothing && error("Filename does not match expected 1D scan pattern: $fname")
    Ei_meV = _fe_parsefloat_token(m.captures[1])
    temperature_K = _fe_parsefloat_token(m.captures[2])
    field_T = _fe_parsefloat_token(m.captures[3])
    qtag = m.captures[4]
    return Ei_meV, temperature_K, field_T, qtag
end

function fe_load_neutron_scan(path::AbstractString)
    Ei_meV, temperature_K, field_T, qtag = fe_parse_neutron_filename(path)
    A = readdlm(path, Float64; comments=true, comment_char='#')
    size(A, 2) >= 6 || error("Expected at least 6 columns in $(basename(path)); got $(size(A, 2))")
    return FEScan1D(
        qtag,
        field_T,
        temperature_K,
        Ei_meV,
        basename(path),
        vec(A[:, 1]),
        vec(A[:, 2]),
        vec(A[:, 3]),
        vec(A[:, 4]),
        vec(A[:, 5]),
        vec(A[:, 6]),
    )
end

function _fe_copy_scan_with_intensity(s::FEScan1D, intensity::Vector{Float64})
    return FEScan1D(
        s.qtag,
        s.field_T,
        s.temperature_K,
        s.Ei_meV,
        s.filename,
        intensity,
        copy(s.error),
        copy(s.energy),
        copy(s.K),
        copy(s.L),
        copy(s.H),
    )
end

function fe_load_neutron_scans(base_dir::AbstractString;
                               Ei_meV::Real,
                               temperature_K::Real,
                               fields_T::Vector{Float64},
                               qtags::Vector{String})
    isdir(base_dir) || error("Could not find neutron 1D scan directory: $base_dir")

    out = Dict{Float64,Dict{String,FEScan1D}}()
    for B in fields_T
        out[Float64(B)] = Dict{String,FEScan1D}()
    end

    files = filter(f -> endswith(f, ".dat") && occursin("_Escan_", basename(f)), readdir(base_dir; join=true))
    for f in files
        s = fe_load_neutron_scan(f)
        abs(s.Ei_meV - Float64(Ei_meV)) <= 1e-6 || continue
        abs(s.temperature_K - Float64(temperature_K)) <= 1e-6 || continue
        s.qtag in qtags || continue
        matches = [B for B in fields_T if abs(B - s.field_T) <= 1e-6]
        isempty(matches) && continue
        out[first(matches)][s.qtag] = s
    end

    for B in fields_T
        for q in qtags
            haskey(out[B], q) || error("Missing neutron scan for field=$(B) T qtag=$q")
        end
    end

    return out
end

function _fe_linear_interpolate(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xnew::AbstractVector{<:Real})
    length(x) == length(y) || error("Interpolation x/y length mismatch")
    p = sortperm(Float64.(x))
    xs = Float64.(x)[p]
    ys = Float64.(y)[p]
    out = similar(Float64.(xnew))
    for (i, xv0) in enumerate(Float64.(xnew))
        xv = Float64(xv0)
        if xv <= xs[1]
            out[i] = ys[1]
        elseif xv >= xs[end]
            out[i] = ys[end]
        else
            j = searchsortedlast(xs, xv)
            t = (xv - xs[j]) / (xs[j + 1] - xs[j])
            out[i] = (1.0 - t) * ys[j] + t * ys[j + 1]
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Neutron background handling
# -----------------------------------------------------------------------------
#
# The default feature-extraction background should match the co-fit logic as
# closely as possible.  The co-fit uses a tail/minimum-field background:
#
#   1. Take the pointwise minimum over 0, 9, and 14 T in a low-energy window
#      and high-energy tail region.
#   2. For K and M cuts, estimate the structured 0 T residual under the broad
#      finite-field peak by fitting a smooth 0 T continuum outside a residual
#      window.
#   3. Add that residual structure back into the smooth background anchor
#      points and interpolate a shape-preserving background over the scan grid.
#
# A simpler :field_minimum option is kept only for diagnostics.

const FE_BG_DEFAULT_LOW_WINDOW = (0.0, 0.75)
const FE_BG_DEFAULT_HIGH_THRESHOLD = 2.5
const FE_BG_DEFAULT_STRUCTURED_FIT_WINDOW = (1.0, 3.0)
const FE_BG_DEFAULT_RESIDUAL_WINDOWS = Dict(
    "0p33_0p33_0" => (1.675, 2.375),
    "0p5_0_0" => (1.825, 2.425),
)

function _fe_common_energy_grid(qdict_by_field::Dict{Float64,FEScan1D}, fields::Vector{Float64}; atol::Real=1e-10)
    for B in fields
        haskey(qdict_by_field, B) || error("Missing required field $(B) T for background construction")
    end
    Eref = qdict_by_field[fields[1]].energy
    for B in fields[2:end]
        E = qdict_by_field[B].energy
        length(E) == length(Eref) || error("Energy grid length mismatch for field $(B) T")
        maximum(abs.(E .- Eref)) <= Float64(atol) || error("Energy grids are not identical for field $(B) T")
    end
    return copy(Eref)
end

function _fe_sort_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || error("x/y length mismatch")
    p = sortperm(Float64.(x))
    return Float64.(x)[p], Float64.(y)[p]
end

struct FEPchipInterpolator
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

function _fe_pchip_endpoint_slope(h1::Float64, h2::Float64, d1::Float64, d2::Float64)
    m = ((2.0*h1 + h2)*d1 - h1*d2) / (h1 + h2)
    if sign(m) != sign(d1)
        return 0.0
    elseif sign(d1) != sign(d2) && abs(m) > abs(3.0*d1)
        return 3.0*d1
    else
        return m
    end
end

function _fe_pchip_interpolator(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real})
    x, y = _fe_sort_xy(xin, yin)
    n = length(x)
    n >= 2 || error("Need at least two points for PCHIP interpolation")
    any(diff(x) .<= 0) && error("PCHIP x-values must be strictly increasing")

    if n == 2
        d = (y[2] - y[1]) / (x[2] - x[1])
        return FEPchipInterpolator(x, y, [d, d])
    end

    h = diff(x)
    d = diff(y) ./ h
    m = zeros(Float64, n)
    m[1] = _fe_pchip_endpoint_slope(h[1], h[2], d[1], d[2])
    m[n] = _fe_pchip_endpoint_slope(h[end], h[end-1], d[end], d[end-1])

    for k in 2:n-1
        if d[k-1] == 0.0 || d[k] == 0.0 || sign(d[k-1]) != sign(d[k])
            m[k] = 0.0
        else
            w1 = 2.0*h[k] + h[k-1]
            w2 = h[k] + 2.0*h[k-1]
            m[k] = (w1 + w2) / (w1/d[k-1] + w2/d[k])
        end
    end

    return FEPchipInterpolator(x, y, m)
end

function (p::FEPchipInterpolator)(x0::Real)
    x = p.x
    y = p.y
    m = p.m
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

(p::FEPchipInterpolator)(xv::AbstractVector{<:Real}) = [p(x) for x in xv]

function _fe_interpolate_background(Egrid::AbstractVector{<:Real}, Eraw::AbstractVector{<:Real}, Iraw::AbstractVector{<:Real}; interpolation_kind::Symbol=:pchip)
    if interpolation_kind == :linear
        return _fe_linear_interpolate(Eraw, Iraw, Egrid)
    elseif interpolation_kind == :pchip
        return _fe_pchip_interpolator(Eraw, Iraw)(Egrid)
    else
        error("Unsupported feature background interpolation_kind=$(interpolation_kind). Use :pchip or :linear.")
    end
end

function _fe_min_over_fields_background(qdict_by_field::Dict{Float64,FEScan1D};
                                        fields::Vector{Float64},
                                        low_window::Tuple{Float64,Float64}=FE_BG_DEFAULT_LOW_WINDOW,
                                        high_threshold::Real=FE_BG_DEFAULT_HIGH_THRESHOLD)
    E = _fe_common_energy_grid(qdict_by_field, fields)
    I_by_field = [qdict_by_field[B].intensity for B in fields]
    Imin = similar(E)
    for i in eachindex(E)
        Imin[i] = minimum(I[i] for I in I_by_field)
    end
    lo, hi = low_window
    raw_mask = ((E .>= lo) .& (E .<= hi)) .| (E .> Float64(high_threshold))
    return E[raw_mask], Imin[raw_mask]
end

function _fe_continuum_design(E::AbstractVector{<:Real}, power::Real; energy_offset::Real=0.15)
    Ef = Float64.(E)
    return hcat(ones(length(Ef)), (Ef .+ Float64(energy_offset)) .^ (-Float64(power)))
end

function _fe_weighted_lsq(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, err::Union{Nothing,AbstractVector{<:Real}})
    Xf = Float64.(X)
    yf = Float64.(y)
    if err !== nothing
        σ = max.(Float64.(err), eps(Float64))
        sw = 1.0 ./ σ
        Xf = Xf .* sw
        yf = yf .* sw
    end
    return (transpose(Xf) * Xf) \ (transpose(Xf) * yf)
end

function _fe_fit_zeroT_power_tail_baseline(Efit::AbstractVector{<:Real}, Ifit::AbstractVector{<:Real}, errfit::Union{Nothing,AbstractVector{<:Real}}=nothing;
                                           power_grid=collect(0.5:0.05:8.0),
                                           energy_offset::Real=0.15,
                                           positive_tail::Bool=true)
    best_score = Inf
    best_coeff = Float64[]
    best_power = NaN
    for power in power_grid
        X = _fe_continuum_design(Efit, power; energy_offset=energy_offset)
        coeff = _fe_weighted_lsq(X, Ifit, errfit)
        if positive_tail && coeff[2] < 0
            continue
        end
        pred = X * coeff
        resid = Float64.(Ifit) .- pred
        score = errfit === nothing ? mean(resid.^2) : mean((resid ./ max.(Float64.(errfit), eps(Float64))).^2)
        if score < best_score
            best_score = score
            best_coeff = Float64.(coeff)
            best_power = Float64(power)
        end
    end
    if isempty(best_coeff)
        return _fe_fit_zeroT_power_tail_baseline(Efit, Ifit, errfit;
            power_grid=power_grid,
            energy_offset=energy_offset,
            positive_tail=false,
        )
    end
    baseline(Enew) = _fe_continuum_design(Enew, best_power; energy_offset=energy_offset) * best_coeff
    return baseline, best_coeff, best_power, best_score
end

function _fe_structured_residual_points(qdict_by_field::Dict{Float64,FEScan1D}, residual_window::Tuple{Float64,Float64};
                                        fit_window::Tuple{Float64,Float64}=FE_BG_DEFAULT_STRUCTURED_FIT_WINDOW,
                                        energy_offset::Real=0.15,
                                        power_grid=collect(0.5:0.05:8.0),
                                        clip_negative_residuals::Bool=false)
    haskey(qdict_by_field, 0.0) || error("Tail background residual modeling requires 0 T data")
    s0 = qdict_by_field[0.0]
    E = s0.energy
    I0 = s0.intensity
    fit_lo, fit_hi = fit_window
    res_lo, res_hi = residual_window
    peak_mask = (E .>= res_lo) .& (E .<= res_hi)
    fit_mask = (E .>= fit_lo) .& (E .<= fit_hi) .& .!peak_mask
    any(fit_mask) || error("No points found for 0 T continuum fit outside residual window $(residual_window)")
    any(peak_mask) || error("No points found inside residual window $(residual_window)")

    baseline_fun, coeff, power, score = _fe_fit_zeroT_power_tail_baseline(
        E[fit_mask], I0[fit_mask], s0.error[fit_mask];
        power_grid=power_grid,
        energy_offset=energy_offset,
    )
    Eres = E[peak_mask]
    residual = I0[peak_mask] .- baseline_fun(Eres)
    if clip_negative_residuals
        residual = max.(residual, 0.0)
    end
    @printf("%s feature background 0T power-tail baseline: power=%.4g score=%.4g coeff=%s\n",
            s0.qtag, power, score, repr(coeff))
    return Eres, residual
end

function _fe_make_tail_background_for_q(qtag::String,
                                        qdict_by_field::Dict{Float64,FEScan1D};
                                        fields::Vector{Float64},
                                        low_window::Tuple{Float64,Float64}=FE_BG_DEFAULT_LOW_WINDOW,
                                        high_threshold::Real=FE_BG_DEFAULT_HIGH_THRESHOLD,
                                        residual_windows::Dict{String,Tuple{Float64,Float64}}=FE_BG_DEFAULT_RESIDUAL_WINDOWS,
                                        structured_fit_window::Tuple{Float64,Float64}=FE_BG_DEFAULT_STRUCTURED_FIT_WINDOW,
                                        interpolation_kind::Symbol=:pchip,
                                        energy_offset::Real=0.15,
                                        clip_negative_residuals::Bool=false)
    Egrid = _fe_common_energy_grid(qdict_by_field, fields)
    E_lowhigh, I_lowhigh = _fe_min_over_fields_background(
        qdict_by_field;
        fields=fields,
        low_window=low_window,
        high_threshold=high_threshold,
    )

    Eraw = copy(E_lowhigh)
    Iraw = copy(I_lowhigh)

    if haskey(residual_windows, qtag)
        Eres, residual = _fe_structured_residual_points(
            qdict_by_field,
            residual_windows[qtag];
            fit_window=structured_fit_window,
            energy_offset=energy_offset,
            clip_negative_residuals=clip_negative_residuals,
        )
        bg_without_residual = _fe_interpolate_background(Egrid, E_lowhigh, I_lowhigh; interpolation_kind=interpolation_kind)
        residual_base = _fe_linear_interpolate(Egrid, bg_without_residual, Eres)
        append!(Eraw, Eres)
        append!(Iraw, residual_base .+ residual)
    end

    bg = _fe_interpolate_background(Egrid, Eraw, Iraw; interpolation_kind=interpolation_kind)
    return Egrid, bg
end

function fe_apply_neutron_background(scans_by_field::Dict{Float64,Dict{String,FEScan1D}};
                                     qtags::Vector{String},
                                     background_mode::Symbol=:tail_bgsub,
                                     background_fields_T::Vector{Float64}=sort(collect(keys(scans_by_field))),
                                     low_window::Tuple{Float64,Float64}=FE_BG_DEFAULT_LOW_WINDOW,
                                     high_threshold::Real=FE_BG_DEFAULT_HIGH_THRESHOLD,
                                     residual_windows::Dict{String,Tuple{Float64,Float64}}=FE_BG_DEFAULT_RESIDUAL_WINDOWS,
                                     structured_fit_window::Tuple{Float64,Float64}=FE_BG_DEFAULT_STRUCTURED_FIT_WINDOW,
                                     interpolation_kind::Symbol=:pchip,
                                     energy_offset::Real=0.15,
                                     clip_negative_residuals::Bool=false)
    if background_mode in (:none, :raw)
        return scans_by_field, Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
    end

    supported = (:tail_bgsub, :spline_bgsub, :cofit_tail, :field_minimum)
    background_mode in supported || error("Unsupported neutron background_mode=$(background_mode). Supported: $(supported), :none, :raw")

    bg = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
    corrected = Dict{Float64,Dict{String,FEScan1D}}()
    for (B, qdict) in scans_by_field
        corrected[B] = Dict{String,FEScan1D}()
    end

    fields = sort(Float64.(background_fields_T))

    for q in qtags
        qdict_by_field = Dict{Float64,FEScan1D}()
        for B in fields
            haskey(scans_by_field, B) && haskey(scans_by_field[B], q) && (qdict_by_field[B] = scans_by_field[B][q])
        end
        isempty(qdict_by_field) && error("No background scans available for qtag=$(q)")

        Eref, Ibg = if background_mode == :field_minimum
            # Diagnostic-only background: pointwise minimum over the full scan.
            available = [qdict_by_field[B] for B in fields if haskey(qdict_by_field, B)]
            E0 = available[1].energy
            Imat = zeros(Float64, length(E0), length(available))
            for (j, s) in enumerate(available)
                Imat[:, j] = length(s.energy) == length(E0) && maximum(abs.(s.energy .- E0)) <= 1e-8 ?
                    s.intensity : _fe_linear_interpolate(s.energy, s.intensity, E0)
            end
            copy(E0), vec(minimum(Imat; dims=2))
        else
            _fe_make_tail_background_for_q(
                q,
                qdict_by_field;
                fields=fields,
                low_window=low_window,
                high_threshold=high_threshold,
                residual_windows=residual_windows,
                structured_fit_window=structured_fit_window,
                interpolation_kind=interpolation_kind,
                energy_offset=energy_offset,
                clip_negative_residuals=clip_negative_residuals,
            )
        end

        bg[q] = (copy(Eref), copy(Ibg))

        for (B, qdict) in scans_by_field
            s = qdict[q]
            bg_on_scan = length(s.energy) == length(Eref) && maximum(abs.(s.energy .- Eref)) <= 1e-8 ?
                Ibg : _fe_linear_interpolate(Eref, Ibg, s.energy)
            corrected[B][q] = _fe_copy_scan_with_intensity(s, s.intensity .- bg_on_scan)
        end
    end

    return corrected, bg
end

# -----------------------------------------------------------------------------
# Neutron feature model
# -----------------------------------------------------------------------------

_fe_gaussian(E, amp, mu, sigma) = amp .* exp.(-0.5 .* ((E .- mu) ./ sigma).^2)
_fe_qsym(qtag::String, stem::String) = Symbol(stem, "_", FE_Q_SHORT[qtag])

function _fe_fit_window_mask(E::AbstractVector{<:Real}, qtag::String, fit_windows_by_q::Dict{String,Vector{Tuple{Float64,Float64}}})
    haskey(fit_windows_by_q, qtag) || error("No neutron fit window for qtag=$qtag")
    mask = falses(length(E))
    for (lo, hi) in fit_windows_by_q[qtag]
        mask .|= (E .>= lo) .& (E .<= hi)
    end
    return mask
end

function fe_neutron_model_for_scan(E::AbstractVector{<:Real}, qtag::String, p::Dict{Symbol,Float64})
    Ed = Float64.(E)
    disp = _fe_gaussian(
        Ed,
        p[_fe_qsym(qtag, "amp_disp")],
        p[_fe_qsym(qtag, "mu_disp")],
        p[_fe_qsym(qtag, "sigma_disp")],
    )
    flat = _fe_gaussian(
        Ed,
        p[_fe_qsym(qtag, "amp_flat")],
        p[:flat_mu_meV],
        p[:flat_sigma_meV],
    )
    offset_value = haskey(p, _fe_qsym(qtag, "offset")) ? p[_fe_qsym(qtag, "offset")] : 0.0
    offset = fill(offset_value, length(Ed))
    total = disp .+ flat .+ offset
    return (; energy_meV=Ed, dispersive=disp, flat=flat, offset=offset, total=total)
end

function _fe_robust_peak_guess(s::FEScan1D, qtag::String, fit_windows_by_q)
    mask = isfinite.(s.energy) .& isfinite.(s.intensity) .& _fe_fit_window_mask(s.energy, qtag, fit_windows_by_q)
    count(mask) >= 3 || return 1e-4, median(s.intensity), 1.5
    E = s.energy[mask]
    y = s.intensity[mask]
    off = median(y)
    imax = argmax(y)
    amp = max(y[imax] - off, 1e-5)
    return amp, off, E[imax]
end

function fe_neutron_feature_specs(data::Dict{String,FEScan1D};
                                  field_T::Real,
                                  qtags::Vector{String},
                                  fit_windows_by_q,
                                  g_flat_guess::Real=3.0,
                                  fit_constant_offset::Bool=true,
                                  cofit_params=nothing,
                                  resolution_sigma_meV::Real=0.08,
                                  amp_split_flat_fraction::Real=0.5)
    B = Float64(field_T)

    use_cofit = cofit_params !== nothing
    flat_mu_guess = if use_cofit
        _fe_flat_center_guess_from_cofit(B, cofit_params)
    else
        Float64(g_flat_guess) * FE_MU_B_MEV_PER_T * B
    end
    flat_mu_guess = clamp(flat_mu_guess, 0.6, 3.6)

    flat_sigma_guess = if use_cofit
        _fe_flat_sigma_guess_from_cofit(B, cofit_params; resolution_sigma_meV=resolution_sigma_meV)
    else
        0.04 + 0.018 * B
    end
    flat_sigma_guess = clamp(flat_sigma_guess, 0.08, 0.95)

    specs = FEParamSpec[
        FEParamSpec(:flat_mu_meV, 0.3, 4.0, flat_mu_guess),
        FEParamSpec(:flat_sigma_meV, 0.03, 1.4, flat_sigma_guess),
    ]

    fflat = clamp(Float64(amp_split_flat_fraction), 0.05, 0.95)

    for q in qtags
        haskey(data, q) || continue
        s = data[q]
        amp0, off0_raw, peakE = _fe_robust_peak_guess(s, q, fit_windows_by_q)
        off0 = fit_constant_offset ? off0_raw : 0.0

        disp_mu_guess = if use_cofit
            _fe_dispersion_guess_from_cofit(q, B, cofit_params)
        else
            peakE
        end
        if !isfinite(disp_mu_guess)
            disp_mu_guess = peakE
        end
        disp_mu_guess = clamp(disp_mu_guess, 0.25, 4.0)

        disp_sigma_guess = if use_cofit
            _fe_dispersion_sigma_guess_from_cofit(q, B, cofit_params; resolution_sigma_meV=resolution_sigma_meV)
        else
            0.20
        end
        disp_sigma_guess = clamp(disp_sigma_guess, 0.05, 0.85)

        # Co-fit-informed amplitude guesses: use local signal near each guessed
        # component energy after the already-applied background subtraction. If
        # the two energies overlap, this is still only an initial condition; the
        # fit itself is free to redistribute amplitude.
        local_disp = max(_fe_nearest_intensity_at(s, disp_mu_guess) - off0, 0.0)
        local_flat = max(_fe_nearest_intensity_at(s, flat_mu_guess) - off0, 0.0)
        if !(isfinite(local_disp) && local_disp > 0)
            local_disp = (1.0 - fflat) * amp0
        end
        if !(isfinite(local_flat) && local_flat > 0)
            local_flat = fflat * amp0
        end

        push!(specs, FEParamSpec(_fe_qsym(q, "amp_disp"), 0.0, 0.05, max(local_disp, 1e-6)))
        push!(specs, FEParamSpec(_fe_qsym(q, "mu_disp"), 0.25, 4.0, disp_mu_guess))
        push!(specs, FEParamSpec(_fe_qsym(q, "sigma_disp"), 0.03, 1.2, disp_sigma_guess))
        push!(specs, FEParamSpec(_fe_qsym(q, "amp_flat"), 0.0, 0.05, max(local_flat, 1e-6)))
        if fit_constant_offset
            push!(specs, FEParamSpec(_fe_qsym(q, "offset"), -0.01, 0.01, off0))
        end
    end

    return specs
end

function fe_fit_neutron_features_for_field(data::Dict{String,FEScan1D};
                                           field_T::Real,
                                           qtags::Vector{String},
                                           fit_windows_by_q,
                                           use_errors::Bool=true,
                                           error_floor::Real=2e-5,
                                           maxiters::Int=3000,
                                           show_trace::Bool=false,
                                           fit_constant_offset::Bool=true,
                                           cofit_params=nothing,
                                           resolution_sigma_meV::Real=0.08,
                                           amp_split_flat_fraction::Real=0.5)
    specs = fe_neutron_feature_specs(
        data;
        field_T=field_T,
        qtags=qtags,
        fit_windows_by_q=fit_windows_by_q,
        fit_constant_offset=fit_constant_offset,
        cofit_params=cofit_params,
        resolution_sigma_meV=resolution_sigma_meV,
        amp_split_flat_fraction=amp_split_flat_fraction,
    )
    u0 = _fe_unconstrained_initial(specs)

    obj = function(u)
        p = _fe_unpack(u, specs)
        ss = 0.0
        n = 0
        for q in qtags
            s = data[q]
            mask = isfinite.(s.energy) .& isfinite.(s.intensity) .& _fe_fit_window_mask(s.energy, q, fit_windows_by_q)
            if use_errors
                mask .&= isfinite.(s.error) .& (s.error .> 0)
            end
            count(mask) >= 3 || continue
            m = fe_neutron_model_for_scan(s.energy[mask], q, p)
            any(!isfinite, m.total) && return 1e30
            denom = use_errors ? sqrt.(s.error[mask].^2 .+ Float64(error_floor)^2) : ones(count(mask))
            r = (s.intensity[mask] .- m.total) ./ denom
            ss += sum(r.^2)
            n += length(r)
        end
        return n > 0 ? ss / n : 1e30
    end

    optres = optimize(obj, u0, NelderMead(), Optim.Options(iterations=maxiters, show_trace=show_trace))
    pbest = _fe_unpack(Optim.minimizer(optres), specs)
    return (; params=pbest, specs=specs, optimizer_result=optres, objective=Optim.minimum(optres))
end

function fe_neutron_feature_rows(field_T::Real, p::Dict{Symbol,Float64}; qtags::Vector{String})
    rows = NamedTuple[]
    Bstr = @sprintf("%.6g", Float64(field_T))
    push!(rows, (block="neutron", field_T=Bstr, qtag="shared", feature="flat_center", value=p[:flat_mu_meV], units="meV", note="shared across q at fixed field"))
    push!(rows, (block="neutron", field_T=Bstr, qtag="shared", feature="flat_sigma", value=p[:flat_sigma_meV], units="meV", note="Gaussian sigma shared across q at fixed field"))
    push!(rows, (block="neutron", field_T=Bstr, qtag="shared", feature="flat_fwhm", value=FE_GAUSS_FWHM * p[:flat_sigma_meV], units="meV", note="2.3548*sigma"))
    if field_T != 0
        push!(rows, (block="neutron", field_T=Bstr, qtag="shared", feature="flat_g_eff_from_center", value=p[:flat_mu_meV] / (FE_MU_B_MEV_PER_T * Float64(field_T)), units="dimensionless", note="E_flat/(mu_B B); empirical feature"))
    end

    for q in qtags
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="dispersive_center", value=p[_fe_qsym(q, "mu_disp")], units="meV", note="independent dispersive Gaussian center"))
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="dispersive_sigma", value=p[_fe_qsym(q, "sigma_disp")], units="meV", note="independent dispersive Gaussian sigma"))
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="dispersive_fwhm", value=FE_GAUSS_FWHM * p[_fe_qsym(q, "sigma_disp")], units="meV", note="2.3548*sigma"))
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="dispersive_height", value=p[_fe_qsym(q, "amp_disp")], units="intensity", note="Gaussian peak height"))
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="flat_height", value=p[_fe_qsym(q, "amp_flat")], units="intensity", note="flat/nondispersive Gaussian peak height"))
        ratio = p[_fe_qsym(q, "amp_disp")] > 0 ? p[_fe_qsym(q, "amp_flat")] / p[_fe_qsym(q, "amp_disp")] : NaN
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="flat_to_dispersive_height_ratio", value=ratio, units="dimensionless", note="empirical intensity ratio"))
        offset_value = haskey(p, _fe_qsym(q, "offset")) ? p[_fe_qsym(q, "offset")] : 0.0
        offset_note = haskey(p, _fe_qsym(q, "offset")) ? "per-cut fitted residual offset" : "fixed to zero by feature-extraction controls"
        push!(rows, (block="neutron", field_T=Bstr, qtag=q, feature="constant_offset", value=offset_value, units="intensity", note=offset_note))
    end
    return rows
end

function fe_neutron_derived_field_rows(neutron_fits::Dict{Float64,NamedTuple}; qtags::Vector{String})
    rows = NamedTuple[]
    fields = sort(collect(keys(neutron_fits)))
    length(fields) < 2 && return rows
    for i in 1:(length(fields)-1)
        B1 = fields[i]
        B2 = fields[i+1]
        p1 = neutron_fits[B1].params
        p2 = neutron_fits[B2].params
        dB = B2 - B1
        dB == 0 && continue
        push!(rows, (block="neutron", field_T=@sprintf("%.6g_to_%.6g", B1, B2), qtag="shared", feature="flat_center_slope_dE_dB", value=(p2[:flat_mu_meV]-p1[:flat_mu_meV])/dB, units="meV/T", note="field derivative from extracted features"))
        push!(rows, (block="neutron", field_T=@sprintf("%.6g_to_%.6g", B1, B2), qtag="shared", feature="flat_g_eff_from_slope", value=(p2[:flat_mu_meV]-p1[:flat_mu_meV])/(FE_MU_B_MEV_PER_T*dB), units="dimensionless", note="dE/dB divided by mu_B"))
        push!(rows, (block="neutron", field_T=@sprintf("%.6g_to_%.6g", B1, B2), qtag="shared", feature="flat_sigma_slope_dsigma_dB", value=(p2[:flat_sigma_meV]-p1[:flat_sigma_meV])/dB, units="meV/T", note="field derivative from extracted features"))
        for q in qtags
            push!(rows, (block="neutron", field_T=@sprintf("%.6g_to_%.6g", B1, B2), qtag=q, feature="dispersive_center_slope_dE_dB", value=(p2[_fe_qsym(q,"mu_disp")]-p1[_fe_qsym(q,"mu_disp")])/dB, units="meV/T", note="field derivative from extracted features"))
        end
    end
    return rows
end

# -----------------------------------------------------------------------------
# Magnetization feature model
# -----------------------------------------------------------------------------

function fe_read_magnetization_csv(filename::AbstractString)
    isfile(filename) || error("Could not find magnetization CSV: $filename")
    B = Float64[]
    M = Float64[]
    open(filename, "r") do io
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
    isempty(B) && error("No numeric rows found in $filename")
    p = sortperm(B)
    return (; B_T=B[p], M_muB_per_Yb=M[p], filename=filename)
end

function fe_filter_magnetization_window(experiment; Bmin_fit_T::Real=0.0, Bmax_fit_T::Real=7.0)
    mask = isfinite.(experiment.B_T) .& isfinite.(experiment.M_muB_per_Yb) .&
           (experiment.B_T .>= Float64(Bmin_fit_T)) .& (experiment.B_T .<= Float64(Bmax_fit_T))
    count(mask) >= 4 || error("Too few magnetization points in fit window")
    return (; B_T=experiment.B_T[mask], M_muB_per_Yb=experiment.M_muB_per_Yb[mask], filename=experiment.filename)
end

function _fe_normal_quadrature(n::Int=101; zmax::Real=5.0)
    n >= 5 || error("normal quadrature needs at least 5 points")
    z = collect(range(-Float64(zmax), Float64(zmax); length=n))
    w = exp.(-0.5 .* z.^2)
    w ./= sum(w)
    return (; z, w)
end

function _fe_spinhalf_tanh(B_fields_T::AbstractVector{<:Real}, g::Real, T_K::Real)
    B = Float64.(B_fields_T)
    T = Float64(T_K)
    T <= 0 && return sign.(B)
    arg = Float64(g) .* FE_MU_B_MEV_PER_T .* B ./ (2.0 * FE_KB_MEV_PER_K * T)
    return tanh.(arg)
end

function _fe_broadened_linear_saturation(B_fields_T::AbstractVector{<:Real}, Bsat_T::Real, sigma_Bsat_T::Real, quad)
    B = Float64.(B_fields_T)
    B0 = Float64(Bsat_T)
    sig = Float64(sigma_Bsat_T)
    if sig <= 1e-8
        return clamp.(B ./ max(B0, eps(Float64)), 0.0, 1.0)
    end

    Bloc_all = B0 .+ sig .* quad.z
    valid = isfinite.(Bloc_all) .& (Bloc_all .> 1e-8)
    count(valid) > 0 || return fill(NaN, length(B))
    Bloc = Bloc_all[valid]
    w = copy(quad.w[valid])
    w ./= sum(w)

    out = zeros(Float64, length(B))
    for (i, b) in enumerate(B)
        ba = max(Float64(b), 0.0)
        total = 0.0
        for j in eachindex(Bloc)
            total += w[j] * clamp(ba / Bloc[j], 0.0, 1.0)
        end
        out[i] = total
    end
    return out
end

function fe_magnetization_model(B_fields_T::AbstractVector{<:Real}, p::Dict{Symbol,Float64}, quad;
                                temperature_K::Real=0.42)
    B = Float64.(B_fields_T)
    offset = fill(p[:M0_offset_muB], length(B))
    vv = p[:chi_vv_muB_per_T] .* B
    flat = p[:M_flat_sat_muB] .* _fe_spinhalf_tanh(B, p[:g_flat], temperature_K)
    disp = p[:M_disp_sat_muB] .* _fe_broadened_linear_saturation(B, p[:Bsat_disp_T], p[:sigma_Bsat_disp_T], quad)
    total = offset .+ vv .+ flat .+ disp
    return (; B_T=B, offset=offset, van_vleck=vv, flat=flat, dispersive=disp, total=total)
end

function fe_magnetization_feature_specs(data;
                                        initial_chi::Real=0.07,
                                        initial_g_flat::Real=3.0)
    Mmax = maximum(data.M_muB_per_Yb)
    return FEParamSpec[
        FEParamSpec(:M0_offset_muB, -0.25, 0.25, 0.0),
        FEParamSpec(:chi_vv_muB_per_T, 0.0, 0.25, Float64(initial_chi)),
        FEParamSpec(:M_flat_sat_muB, 0.0, 1.50, 0.15 * max(Mmax, 0.1)),
        FEParamSpec(:g_flat, 0.5, 8.0, Float64(initial_g_flat)),
        FEParamSpec(:M_disp_sat_muB, 0.0, 2.50, 0.85 * max(Mmax, 0.1)),
        FEParamSpec(:Bsat_disp_T, 0.25, 20.0, 4.0),
        FEParamSpec(:sigma_Bsat_disp_T, 0.0, 10.0, 1.0),
    ]
end

function fe_fit_magnetization_features(data;
                                       temperature_K::Real=0.42,
                                       initial_chi::Real=0.07,
                                       initial_g_flat::Real=3.0,
                                       quad_n::Int=101,
                                       maxiters::Int=3000,
                                       show_trace::Bool=false)
    quad = _fe_normal_quadrature(quad_n)
    specs = fe_magnetization_feature_specs(data; initial_chi=initial_chi, initial_g_flat=initial_g_flat)
    u0 = _fe_unconstrained_initial(specs)
    sigma = max(0.01 * (maximum(data.M_muB_per_Yb) - minimum(data.M_muB_per_Yb)), 1e-4)

    obj = function(u)
        p = _fe_unpack(u, specs)
        m = fe_magnetization_model(data.B_T, p, quad; temperature_K=temperature_K)
        any(!isfinite, m.total) && return 1e30
        r = (data.M_muB_per_Yb .- m.total) ./ sigma
        return mean(r.^2)
    end

    optres = optimize(obj, u0, NelderMead(), Optim.Options(iterations=maxiters, show_trace=show_trace))
    pbest = _fe_unpack(Optim.minimizer(optres), specs)
    return (; params=pbest, specs=specs, optimizer_result=optres, objective=Optim.minimum(optres), quad=quad)
end

function fe_magnetization_feature_rows(p::Dict{Symbol,Float64})
    rows = NamedTuple[]
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="M0_offset", value=p[:M0_offset_muB], units="mu_B/Yb", note="constant offset"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="chi_vv", value=p[:chi_vv_muB_per_T], units="mu_B/Yb/T", note="linear Van Vleck-like slope"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="M_flat_sat", value=p[:M_flat_sat_muB], units="mu_B/Yb", note="flat/nondispersive tanh saturation amplitude"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="g_flat", value=p[:g_flat], units="dimensionless", note="flat/nondispersive empirical tanh g-factor"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="M_disp_sat", value=p[:M_disp_sat_muB], units="mu_B/Yb", note="dispersive broadened-saturation amplitude"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="Bsat_disp", value=p[:Bsat_disp_T], units="T", note="mean saturation field for dispersive empirical component"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="sigma_Bsat_disp", value=p[:sigma_Bsat_disp_T], units="T", note="Gaussian broadening of local saturation field"))
    push!(rows, (block="magnetization", field_T="all", qtag="all", feature="dispersive_initial_slope", value=p[:M_disp_sat_muB] / max(p[:Bsat_disp_T], eps(Float64)), units="mu_B/Yb/T", note="M_disp_sat/Bsat_disp"))
    return rows
end

# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------

function _fe_csv_cell(x)
    s = string(x)
    s = replace(s, '"' => "\"\"")
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * s * "\""
    end
    return s
end

function fe_write_feature_summary_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "block,field_T,qtag,feature,value,units,note")
        for r in rows
            println(io, join((_fe_csv_cell(r.block), _fe_csv_cell(r.field_T), _fe_csv_cell(r.qtag), _fe_csv_cell(r.feature), @sprintf("%.12g", r.value), _fe_csv_cell(r.units), _fe_csv_cell(r.note)), ","))
        end
    end
end

function fe_write_neutron_components_csv(path::AbstractString,
                                         scans_by_field::Dict{Float64,Dict{String,FEScan1D}},
                                         neutron_fits::Dict{Float64,NamedTuple};
                                         qtags::Vector{String},
                                         fit_windows_by_q)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "field_T,qtag,energy_meV,data_intensity,data_error,dispersive,flat,offset,total,in_fit_window")
        for B in sort(collect(keys(neutron_fits)))
            p = neutron_fits[B].params
            for q in qtags
                s = scans_by_field[B][q]
                m = fe_neutron_model_for_scan(s.energy, q, p)
                fitmask = _fe_fit_window_mask(s.energy, q, fit_windows_by_q)
                for i in eachindex(s.energy)
                    println(io, join((@sprintf("%.12g", B), q, @sprintf("%.12g", s.energy[i]), @sprintf("%.12g", s.intensity[i]), @sprintf("%.12g", s.error[i]), @sprintf("%.12g", m.dispersive[i]), @sprintf("%.12g", m.flat[i]), @sprintf("%.12g", m.offset[i]), @sprintf("%.12g", m.total[i]), fitmask[i]), ","))
                end
            end
        end
    end
end

function fe_write_neutron_background_csv(path::AbstractString, bg)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "qtag,energy_meV,background_intensity")
        for q in sort(collect(keys(bg)))
            E, I = bg[q]
            for i in eachindex(E)
                println(io, join((q, @sprintf("%.12g", E[i]), @sprintf("%.12g", I[i])), ","))
            end
        end
    end
end

function fe_write_magnetization_components_csv(path::AbstractString,
                                               experiment_all,
                                               mag_fit;
                                               temperature_K::Real=0.42,
                                               Bmin_plot_T::Real=0.0,
                                               Bmax_plot_T::Real=7.0,
                                               dB_plot_T::Real=0.01)
    mkpath(dirname(path))
    Bgrid = collect(Float64(Bmin_plot_T):Float64(dB_plot_T):Float64(Bmax_plot_T))
    m = fe_magnetization_model(Bgrid, mag_fit.params, mag_fit.quad; temperature_K=temperature_K)
    open(path, "w") do io
        println(io, "B_T,M_total_model,M_disp,M_flat,M_van_vleck,M_offset")
        for i in eachindex(Bgrid)
            println(io, join((@sprintf("%.12g", m.B_T[i]), @sprintf("%.12g", m.total[i]), @sprintf("%.12g", m.dispersive[i]), @sprintf("%.12g", m.flat[i]), @sprintf("%.12g", m.van_vleck[i]), @sprintf("%.12g", m.offset[i])), ","))
        end
    end
end

function fe_plot_feature_overview(path::AbstractString,
                                  scans_by_field::Dict{Float64,Dict{String,FEScan1D}},
                                  neutron_fits::Dict{Float64,NamedTuple},
                                  qtags::Vector{String},
                                  mag_data,
                                  mag_fit;
                                  fit_windows_by_q,
                                  magnetization_temperature_K::Real=0.42,
                                  neutron_ylim::Union{Nothing,Tuple{Float64,Float64}}=nothing)
    fields = sort(collect(keys(neutron_fits)))
    nrows = length(fields) + 1
    ncols = length(qtags)
    fig = Figure(size=(420 * ncols, 300 * nrows))

    for (iB, B) in enumerate(fields)
        p = neutron_fits[B].params
        for (iq, q) in enumerate(qtags)
            s = scans_by_field[B][q]
            ax = Axis(fig[iB, iq], xlabel="Energy transfer (meV)", ylabel="Intensity", title="$(B) T, $(get(FE_Q_LABEL, q, q))")
            if neutron_ylim !== nothing
                ylims!(ax, neutron_ylim[1], neutron_ylim[2])
            end
            scatter!(ax, s.energy, s.intensity; markersize=5)
            m = fe_neutron_model_for_scan(s.energy, q, p)
            lines!(ax, s.energy, m.total; linewidth=2)
            lines!(ax, s.energy, m.dispersive; linewidth=1, linestyle=:dash)
            lines!(ax, s.energy, m.flat; linewidth=1, linestyle=:dot)
            for (lo, hi) in fit_windows_by_q[q]
                vspan!(ax, lo, hi; alpha=0.08)
            end
        end
    end

    axm = Axis(fig[nrows, 1:ncols], xlabel="Field (T)", ylabel="M (μB/Yb)", title="Magnetization empirical feature model")
    scatter!(axm, mag_data.B_T, mag_data.M_muB_per_Yb; markersize=5)
    Bgrid = collect(range(minimum(mag_data.B_T), maximum(mag_data.B_T); length=400))
    mm = fe_magnetization_model(Bgrid, mag_fit.params, mag_fit.quad; temperature_K=magnetization_temperature_K)
    lines!(axm, Bgrid, mm.total; linewidth=2)
    lines!(axm, Bgrid, mm.dispersive; linewidth=1, linestyle=:dash)
    lines!(axm, Bgrid, mm.flat; linewidth=1, linestyle=:dot)
    lines!(axm, Bgrid, mm.van_vleck .+ mm.offset; linewidth=1, linestyle=:dashdot)

    save(path, fig)
    return fig
end

# -----------------------------------------------------------------------------
# Config parsing and top-level workflow
# -----------------------------------------------------------------------------

function _fe_toml_symbol(value)
    value isa Symbol && return value
    value isa AbstractString && return Symbol(value)
    error("Cannot convert $value to Symbol")
end

function _fe_string_vector(v)
    return String[string(x) for x in v]
end

function _fe_float_vector(v)
    return Float64[Float64(x) for x in v]
end

function _fe_tuple2(v; default=nothing)
    if v === nothing
        default === nothing && error("Expected a two-element value, got nothing")
        return default
    end
    length(v) == 2 || error("Expected two-element value, got $(v)")
    return (Float64(v[1]), Float64(v[2]))
end

function _fe_optional_tuple2(v)
    v === nothing && return nothing
    return _fe_tuple2(v)
end

function _fe_residual_windows_from_controls(ctrl)
    if !haskey(ctrl["neutron"], "background") || !haskey(ctrl["neutron"]["background"], "residual_windows")
        return FE_BG_DEFAULT_RESIDUAL_WINDOWS
    end
    table = ctrl["neutron"]["background"]["residual_windows"]
    out = Dict{String,Tuple{Float64,Float64}}()
    for (q, w) in table
        out[string(q)] = _fe_tuple2(w)
    end
    return out
end

function _fe_fit_windows_from_controls(ctrl, qtags::Vector{String})
    table = ctrl["neutron"]["fit_windows"]
    out = Dict{String,Vector{Tuple{Float64,Float64}}}()
    for q in qtags
        haskey(table, q) || error("Missing [neutron.fit_windows] entry for qtag=$q")
        out[q] = Tuple{Float64,Float64}[(Float64(w[1]), Float64(w[2])) for w in table[q]]
    end
    return out
end


function _fe_neutron_initialization_controls(ctrl)
    return get(ctrl["neutron"], "initialization", Dict{String,Any}())
end

function _fe_load_cofit_initialization_params(repo_root::AbstractString, ctrl)
    init_ctrl = _fe_neutron_initialization_controls(ctrl)
    mode = _fe_toml_symbol(get(init_ctrl, "mode", "robust_data"))
    if mode in (:none, :robust_data, :data)
        return nothing
    elseif mode == :cofit_informed
        relpath = String(get(init_ctrl, "best_fit_parameters_toml", "configs/best_fit_parameters.toml"))
        path = normpath(joinpath(repo_root, relpath))
        println("Loading co-fit-informed neutron feature initial guesses from:")
        println(path)
        return load_canonical_model_parameters(path)
    else
        error("Unsupported [neutron.initialization] mode=$(mode). Use robust_data or cofit_informed.")
    end
end

function run_neutron_magnetization_feature_extraction(; repo_root::AbstractString, controls::Dict)
    neutron_dir = normpath(joinpath(repo_root, controls["data"]["neutron_1d_subdir"]))
    magnetization_csv = normpath(joinpath(repo_root, controls["data"]["magnetization_csv"]))
    out_table_dir = normpath(joinpath(repo_root, controls["output"]["feature_table_subdir"]))
    out_figure_dir = normpath(joinpath(repo_root, controls["output"]["figure_subdir"]))
    mkpath(out_table_dir)
    mkpath(out_figure_dir)

    qtags = _fe_string_vector(controls["neutron"]["qtags"])
    fields_T = _fe_float_vector(controls["neutron"]["fields_T"])
    background_fields_T = haskey(controls["neutron"], "background_fields_T") ? _fe_float_vector(controls["neutron"]["background_fields_T"]) : fields_T
    loaded_fields_T = sort(unique(vcat(fields_T, background_fields_T)))
    fit_windows_by_q = _fe_fit_windows_from_controls(controls, qtags)
    init_controls = _fe_neutron_initialization_controls(controls)
    cofit_init_params = _fe_load_cofit_initialization_params(repo_root, controls)
    fit_constant_offset = Bool(get(controls["neutron"], "fit_constant_offset", true))
    resolution_sigma_meV = Float64(get(init_controls, "resolution_sigma_meV", 0.08))
    amp_split_flat_fraction = Float64(get(init_controls, "amp_split_flat_fraction", 0.5))

    println("Neutron feature extraction options:")
    println("  fit_constant_offset = ", fit_constant_offset)
    println("  initialization mode = ", get(init_controls, "mode", "robust_data"))
    println("  resolution_sigma_meV for initial guesses = ", resolution_sigma_meV)
    println()

    println("Loading neutron scans from:")
    println(neutron_dir)
    scans_raw = fe_load_neutron_scans(
        neutron_dir;
        Ei_meV=Float64(controls["neutron"]["Ei_meV"]),
        temperature_K=Float64(controls["neutron"]["temperature_K"]),
        fields_T=loaded_fields_T,
        qtags=qtags,
    )

    background_mode = _fe_toml_symbol(controls["neutron"]["background_mode"])
    bg_controls = get(controls["neutron"], "background", Dict{String,Any}())
    scans, bg = fe_apply_neutron_background(
        scans_raw;
        qtags=qtags,
        background_mode=background_mode,
        background_fields_T=background_fields_T,
        low_window=_fe_tuple2(get(bg_controls, "low_window", [FE_BG_DEFAULT_LOW_WINDOW[1], FE_BG_DEFAULT_LOW_WINDOW[2]])),
        high_threshold=Float64(get(bg_controls, "high_threshold", FE_BG_DEFAULT_HIGH_THRESHOLD)),
        residual_windows=_fe_residual_windows_from_controls(controls),
        structured_fit_window=_fe_tuple2(get(bg_controls, "structured_fit_window", [FE_BG_DEFAULT_STRUCTURED_FIT_WINDOW[1], FE_BG_DEFAULT_STRUCTURED_FIT_WINDOW[2]])),
        interpolation_kind=_fe_toml_symbol(get(bg_controls, "interpolation_kind", "pchip")),
        energy_offset=Float64(get(bg_controls, "energy_offset", 0.15)),
        clip_negative_residuals=Bool(get(bg_controls, "clip_negative_residuals", false)),
    )

    neutron_fits = Dict{Float64,NamedTuple}()
    for B in fields_T
        println("Fitting neutron empirical features for field $(B) T")
        neutron_fits[B] = fe_fit_neutron_features_for_field(
            scans[B];
            field_T=B,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            use_errors=Bool(controls["neutron"]["use_errors"]),
            error_floor=Float64(controls["neutron"]["error_floor"]),
            maxiters=Int(controls["neutron"]["maxiters"]),
            show_trace=Bool(get(controls["neutron"], "show_trace", false)),
            fit_constant_offset=fit_constant_offset,
            cofit_params=cofit_init_params,
            resolution_sigma_meV=resolution_sigma_meV,
            amp_split_flat_fraction=amp_split_flat_fraction,
        )
        println("  objective = ", neutron_fits[B].objective)
    end

    println("Loading magnetization data from:")
    println(magnetization_csv)
    mag_all = fe_read_magnetization_csv(magnetization_csv)
    mag_fit_data = fe_filter_magnetization_window(
        mag_all;
        Bmin_fit_T=Float64(controls["magnetization"]["Bmin_fit_T"]),
        Bmax_fit_T=Float64(controls["magnetization"]["Bmax_fit_T"]),
    )
    println("Fitting magnetization empirical features")
    mag_fit = fe_fit_magnetization_features(
        mag_fit_data;
        temperature_K=Float64(controls["magnetization"]["temperature_K"]),
        initial_chi=Float64(controls["magnetization"]["initial_chi_vv_muB_per_T"]),
        initial_g_flat=Float64(controls["magnetization"]["initial_g_flat"]),
        quad_n=Int(controls["magnetization"]["quad_n"]),
        maxiters=Int(controls["magnetization"]["maxiters"]),
        show_trace=Bool(get(controls["magnetization"], "show_trace", false)),
    )
    println("  objective = ", mag_fit.objective)

    rows = NamedTuple[]
    for B in fields_T
        append!(rows, fe_neutron_feature_rows(B, neutron_fits[B].params; qtags=qtags))
    end
    append!(rows, fe_neutron_derived_field_rows(neutron_fits; qtags=qtags))
    append!(rows, fe_magnetization_feature_rows(mag_fit.params))

    summary_path = joinpath(out_table_dir, "feature_summary.csv")
    neutron_components_path = joinpath(out_table_dir, "neutron_components.csv")
    neutron_background_path = joinpath(out_table_dir, "neutron_background.csv")
    magnetization_components_path = joinpath(out_table_dir, "magnetization_components.csv")

    fe_write_feature_summary_csv(summary_path, rows)
    fe_write_neutron_components_csv(neutron_components_path, scans, neutron_fits; qtags=qtags, fit_windows_by_q=fit_windows_by_q)
    fe_write_neutron_background_csv(neutron_background_path, bg)
    fe_write_magnetization_components_csv(
        magnetization_components_path,
        mag_all,
        mag_fit;
        temperature_K=Float64(controls["magnetization"]["temperature_K"]),
        Bmin_plot_T=Float64(controls["magnetization"]["Bmin_plot_T"]),
        Bmax_plot_T=Float64(controls["magnetization"]["Bmax_plot_T"]),
        dB_plot_T=Float64(controls["magnetization"]["dB_plot_T"]),
    )

    fig_path = nothing
    if Bool(controls["plotting"]["make_plots"])
        fig_path = joinpath(out_figure_dir, "neutron_magnetization_feature_extraction.png")
        fe_plot_feature_overview(
            fig_path,
            scans,
            neutron_fits,
            qtags,
            mag_fit_data,
            mag_fit;
            fit_windows_by_q=fit_windows_by_q,
            magnetization_temperature_K=Float64(controls["magnetization"]["temperature_K"]),
            neutron_ylim=_fe_optional_tuple2(get(controls["plotting"], "neutron_ylim", nothing)),
        )
    end

    println()
    println("Feature extraction complete.")
    println("Feature summary: ", summary_path)
    println("Neutron components: ", neutron_components_path)
    println("Magnetization components: ", magnetization_components_path)
    if fig_path !== nothing
        println("Overview figure: ", fig_path)
    end

    return (; rows, neutron_fits, magnetization_fit=mag_fit, summary_path, neutron_components_path, magnetization_components_path, figure_path=fig_path)
end
