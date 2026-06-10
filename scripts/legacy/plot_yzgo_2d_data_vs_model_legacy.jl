# =============================================================================
# YZGO standalone co-fit: neutron 1D scans + full-field magnetization, shared non-dispersive fraction
# =============================================================================
#
# This standalone file embeds the two original source scripts in private modules
# and adds a co-fit driver at the bottom with one shared non-dispersive relative weight. It does not auto-run the original
# individual fits. Running this file directly runs the co-fit.
#
# Default use:
#   julia YZGO_cofit_neutron_magnetization_two_phase_standalone.jl
#
# Interactive use:
#   include("YZGO_cofit_neutron_magnetization_two_phase_standalone.jl")
#   out = run_yzgo_neutron_magnetization_cofit_shared_fraction(neutron_weight=1.0, magnetization_weight=1.0)
#
# =============================================================================

module YZGONeutronFit
# =============================================================================
# YZGO monolithic 1D fit: two-kernel model
#
# This version removes the Q=(0,1,0)-specific relative scale factor.
# Instead it fits a second, non-dispersive scattering kernel with:
#   - effective S = 1/2
#   - J1 = J2 = 0
#   - fitted gzz2 and sigma_gzz2
#   - fitted intensity relative to the dispersive single-magnon kernel
#
# Plots show: data, dispersive kernel, non-dispersive kernel, and total.
# =============================================================================



# =============================================================================
# Section copied from yzgo_plot_1d_scans.jl, with auto-run block removed
# =============================================================================

# yzgo_plot_1d_scans.jl
# Load, plot, and background-subtract YZGO 1D energy scans at several field values
# and Q centerings.
#
# Usage from Julia:
#   include("yzgo_plot_1d_scans.jl")
#
#   scans = load_scans(BASEDIR)
#   print_scan_summary(scans)
#
#   # Final requested figure: raw top row, model-background-subtracted bottom row.
#   fig, scans_bgsub, background_models = plot_raw_and_tail_bgsubtracted(scans)
#
# The loaded data remain in memory as `scans`, and the corrected data remain in memory
# as `scans_bgsub`, with access like:
#   s = scans["0_1_0"][9.0]
#   s.energy, s.intensity, s.error
#
#   c = scans_bgsub["0_1_0"][9.0]
#   c.energy, c.intensity, c.error
#
#   bg = background_models["0p33_0p33_0"]
#   bg.energy_grid, bg.background

using DelimitedFiles
using Printf
using Statistics
using LinearAlgebra
using GLMakie

const BASEDIR = raw"C:\Users\vdp\ORNL Dropbox\Daniel Pajerowski\YZGO\CNCS_data\SYM_1d_scans"
const OUTDIR = joinpath(BASEDIR, "plots")

struct HistMeta
    label::String
    de_range::Tuple{Float64,Float64}
    de_step::Float64
    h_range::Tuple{Float64,Float64}
    k_range::Tuple{Float64,Float64}
    l_range::Tuple{Float64,Float64}
    plot_ylim::Tuple{Float64,Float64}
end

const Q_ORDER = ["0_1_0", "0p33_0p33_0", "0p5_0_0"]
const FIELD_ORDER = [0.0, 9.0, 14.0]

const QMETA = Dict(
    "0_1_0" => HistMeta(
        "Q = (0, 1, 0)",
        (0.0, 3.3),
        0.05,
        (-0.1, 0.1),
        (0.9, 1.1),
        (-0.3, 0.3),
        (0.0, 0.004),
    ),
    "0p33_0p33_0" => HistMeta(
        "Q = (0.33, 0.33, 0)",
        (0.0, 4.0),
        0.05,
        (0.23, 0.43),
        (0.23, 0.43),
        (-0.3, 0.3),
        (0.0, 0.004),
    ),
    "0p5_0_0" => HistMeta(
        "Q = (0.5, 0, 0)",
        (0.0, 4.0),
        0.05,
        (0.4, 0.6),
        (-0.1, 0.1),
        (-0.3, 0.3),
        (0.0, 0.004),
    ),
)

struct Scan1D
    qtag::String
    field_T::Float64
    temperature_K::Float64
    Ei_meV::Float64
    filename::String
    meta::HistMeta
    intensity::Vector{Float64}
    error::Vector{Float64}
    energy::Vector{Float64}
    K::Vector{Float64}
    L::Vector{Float64}
    H::Vector{Float64}
end

struct NaturalCubicSpline
    x::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

struct BackgroundModel
    qtag::String
    energy_raw::Vector{Float64}
    intensity_raw::Vector{Float64}
    intensity_raw_smoothed::Vector{Float64}
    energy_grid::Vector{Float64}
    background::Vector{Float64}
    lowhigh_energy::Vector{Float64}
    lowhigh_intensity::Vector{Float64}
    residual_window::Union{Nothing,Tuple{Float64,Float64}}
    residual_energy::Vector{Float64}
    residual::Vector{Float64}
    residual_absolute_background::Vector{Float64}
    zeroT_fit_energy::Vector{Float64}
    zeroT_fit_intensity::Vector{Float64}
    zeroT_baseline_in_residual_window::Vector{Float64}
end

parsefloat_token(s::AbstractString) = parse(Float64, replace(s, "p" => "."))

function parse_filename(path::AbstractString)
    fname = basename(path)
    m = match(r"^yzgo_(\d+p\d+)meV_(\d+(?:p\d+)?)K_(\d+(?:p\d+)?)T_Escan_(.+)_SYM\.dat$", fname)
    m === nothing && error("Filename does not match expected pattern: $fname")
    Ei_meV = parsefloat_token(m.captures[1])
    temperature_K = parsefloat_token(m.captures[2])
    field_T = parsefloat_token(m.captures[3])
    qtag = m.captures[4]
    haskey(QMETA, qtag) || error("No histogram metadata defined for qtag = $qtag")
    return Ei_meV, temperature_K, field_T, qtag
end

function load_scan(path::AbstractString)
    Ei_meV, temperature_K, field_T, qtag = parse_filename(path)

    # File columns:
    #   Intensity, Error, DeltaE, [0,K,0], [0,0,L], [H,0,0]
    A = readdlm(path, Float64; comments=true, comment_char='#')
    size(A, 2) >= 6 || error("Expected at least 6 columns in $(basename(path)); got $(size(A, 2))")

    return Scan1D(
        qtag,
        field_T,
        temperature_K,
        Ei_meV,
        basename(path),
        QMETA[qtag],
        vec(A[:, 1]),
        vec(A[:, 2]),
        vec(A[:, 3]),
        vec(A[:, 4]),
        vec(A[:, 5]),
        vec(A[:, 6]),
    )
end

function load_scans(base_dir::AbstractString=BASEDIR)
    files = filter(readdir(base_dir; join=true)) do f
        endswith(f, ".dat") && occursin("_Escan_", basename(f)) && occursin("_SYM.dat", basename(f))
    end

    scans = Dict{String,Dict{Float64,Scan1D}}()
    for f in sort(files)
        s = load_scan(f)
        byfield = get!(scans, s.qtag, Dict{Float64,Scan1D}())
        byfield[s.field_T] = s
    end
    return scans
end

function print_scan_summary(scans::Dict{String,Dict{Float64,Scan1D}})
    for qtag in Q_ORDER
        haskey(scans, qtag) || continue
        println("\n", QMETA[qtag].label)
        for B in sort(collect(keys(scans[qtag])))
            s = scans[qtag][B]
            @printf("  %4.1f T: %3d points, ΔE = %.3f to %.3f meV, I = %.4g to %.4g\n",
                    B, length(s.energy), minimum(s.energy), maximum(s.energy), minimum(s.intensity), maximum(s.intensity))
        end
    end
end

function plot_initial_scans(scans::Dict{String,Dict{Float64,Scan1D}};
                            outdir::AbstractString=OUTDIR,
                            save_png::Bool=true)
    mkpath(outdir)

    fig = Figure(size=(1450, 430))

    for (icol, qtag) in enumerate(Q_ORDER)
        haskey(scans, qtag) || continue
        meta = QMETA[qtag]
        ax = Axis(fig[1, icol];
            title=meta.label,
            xlabel="ΔE (meV)",
            ylabel=icol == 1 ? "Intensity (arb. units)" : "",
            limits=(meta.de_range[1], meta.de_range[2], meta.plot_ylim[1], meta.plot_ylim[2]),
        )

        for B in sort(collect(keys(scans[qtag])))
            s = scans[qtag][B]
            lines!(ax, s.energy, s.intensity; label=@sprintf("%.0f T", B), linewidth=2)
            errorbars!(ax, s.energy, s.intensity, s.error; whiskerwidth=4)
        end

        axislegend(ax; position=:rt, framevisible=false)
    end

    if save_png
        save(joinpath(outdir, "YZGO_1d_Escan_initial_field_overlay.png"), fig)
    end
    return fig
end

# -------------------------------------------------------------------------
# Spline/minimum-field background subtraction
# -------------------------------------------------------------------------

const MIN_BG_LOW_WINDOW = (0.0, 0.75)      # use min over 0, 9, 14 T for 0 <= E <= 0.75 meV
const MIN_BG_HIGH_THRESHOLD = 2.5          # use min over 0, 9, 14 T for E > 2.5 meV
const STRUCTURED_FIT_WINDOW = (1.0, 3.0)   # fit 0 T data in this broad window, excluding peak window

const STRUCTURED_RESIDUAL_WINDOWS = Dict(
    "0p33_0p33_0" => (1.675, 2.375),
    "0p5_0_0" => (1.825, 2.425),
)

function energy_window_mask(E::AbstractVector{<:Real}, windows::Vector{Tuple{Float64,Float64}})
    mask = falses(length(E))
    for (lo, hi) in windows
        mask .|= (E .>= lo) .& (E .<= hi)
    end
    return mask
end

function assert_required_fields(byfield::Dict{Float64,Scan1D}; fields=FIELD_ORDER)
    missing_fields = [B for B in fields if !haskey(byfield, B)]
    isempty(missing_fields) || error("Missing required field scans: $(missing_fields)")
    return nothing
end

function common_energy_grid(byfield::Dict{Float64,Scan1D}; fields=FIELD_ORDER, atol=1e-10)
    assert_required_fields(byfield; fields=fields)
    Eref = byfield[fields[1]].energy
    for B in fields[2:end]
        E = byfield[B].energy
        length(E) == length(Eref) || error("Energy grid length mismatch for $(byfield[B].qtag), field $(B) T")
        maximum(abs.(E .- Eref)) <= atol || error("Energy grids are not identical for $(byfield[B].qtag)")
    end
    return copy(Eref)
end

function min_over_fields_background_raw(byfield::Dict{Float64,Scan1D};
                                        fields=FIELD_ORDER,
                                        low_window=MIN_BG_LOW_WINDOW,
                                        high_threshold=MIN_BG_HIGH_THRESHOLD)
    E = common_energy_grid(byfield; fields=fields)
    I_by_field = [byfield[B].intensity for B in fields]

    Imin = similar(E)
    for i in eachindex(E)
        Imin[i] = minimum(I[i] for I in I_by_field)
    end

    lo, hi = low_window
    raw_mask = ((E .>= lo) .& (E .<= hi)) .| (E .> high_threshold)
    return E[raw_mask], Imin[raw_mask]
end

function sort_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || error("x and y must have the same length")
    p = sortperm(x)
    return Float64.(x[p]), Float64.(y[p])
end

function gaussian_smooth_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, sigma::Real)
    length(x) == length(y) || error("x and y must have same length")
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

function natural_cubic_spline(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real})
    x, y = sort_xy(xin, yin)
    n = length(x)
    n == length(y) || error("x and y must have same length")
    n >= 2 || error("Need at least two points for spline")
    any(diff(x) .<= 0) && error("Spline x-values must be strictly increasing")

    h = diff(x)
    α = zeros(Float64, n)
    for i in 2:n-1
        α[i] = 3.0 / h[i] * (y[i+1] - y[i]) - 3.0 / h[i-1] * (y[i] - y[i-1])
    end

    l = ones(Float64, n)
    μ = zeros(Float64, n)
    z = zeros(Float64, n)
    for i in 2:n-1
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end

    a = y[1:n-1]
    b = zeros(Float64, n-1)
    c = zeros(Float64, n)
    d = zeros(Float64, n-1)
    for j in (n-1):-1:1
        c[j] = z[j] - μ[j] * c[j+1]
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2.0 * c[j]) / 3.0
        d[j] = (c[j+1] - c[j]) / (3.0 * h[j])
    end

    return NaturalCubicSpline(x, a, b, c[1:n-1], d)
end

function (s::NaturalCubicSpline)(x0::Real)
    x = s.x
    n = length(x)
    j = searchsortedlast(x, Float64(x0))
    j = clamp(j, 1, n - 1)
    dx = Float64(x0) - x[j]
    return s.a[j] + s.b[j] * dx + s.c[j] * dx^2 + s.d[j] * dx^3
end

(s::NaturalCubicSpline)(xv::AbstractVector{<:Real}) = [s(x) for x in xv]


struct PchipInterpolator
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

function pchip_endpoint_slope(h1::Float64, h2::Float64, d1::Float64, d2::Float64)
    m = ((2.0*h1 + h2)*d1 - h1*d2) / (h1 + h2)
    if sign(m) != sign(d1)
        return 0.0
    elseif sign(d1) != sign(d2) && abs(m) > abs(3.0*d1)
        return 3.0*d1
    else
        return m
    end
end

function pchip_interpolator(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real})
    x, y = sort_xy(xin, yin)
    n = length(x)
    n == length(y) || error("x and y must have same length")
    n >= 2 || error("Need at least two points for interpolation")
    any(diff(x) .<= 0) && error("Interpolation x-values must be strictly increasing")

    if n == 2
        d = (y[2] - y[1]) / (x[2] - x[1])
        return PchipInterpolator(x, y, [d, d])
    end

    h = diff(x)
    d = diff(y) ./ h
    m = zeros(Float64, n)

    m[1] = pchip_endpoint_slope(h[1], h[2], d[1], d[2])
    m[n] = pchip_endpoint_slope(h[end], h[end-1], d[end], d[end-1])

    for k in 2:n-1
        if d[k-1] == 0.0 || d[k] == 0.0 || sign(d[k-1]) != sign(d[k])
            m[k] = 0.0
        else
            w1 = 2.0*h[k] + h[k-1]
            w2 = h[k] + 2.0*h[k-1]
            m[k] = (w1 + w2) / (w1/d[k-1] + w2/d[k])
        end
    end

    return PchipInterpolator(x, y, m)
end

function (p::PchipInterpolator)(x0::Real)
    x = p.x
    y = p.y
    m = p.m
    n = length(x)
    j = searchsortedlast(x, Float64(x0))
    j = clamp(j, 1, n - 1)
    h = x[j+1] - x[j]
    t = (Float64(x0) - x[j]) / h

    h00 = (2.0*t^3 - 3.0*t^2 + 1.0)
    h10 = (t^3 - 2.0*t^2 + t)
    h01 = (-2.0*t^3 + 3.0*t^2)
    h11 = (t^3 - t^2)

    return h00*y[j] + h10*h*m[j] + h01*y[j+1] + h11*h*m[j+1]
end

(p::PchipInterpolator)(xv::AbstractVector{<:Real}) = [p(x) for x in xv]

function linear_interpolate(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real}, xgrid::AbstractVector{<:Real})
    x, y = sort_xy(xin, yin)
    n = length(x)
    n >= 2 || error("Need at least two points for linear interpolation")
    any(diff(x) .<= 0) && error("Interpolation x-values must be strictly increasing")

    out = Vector{Float64}(undef, length(xgrid))
    for (i, x0r) in enumerate(xgrid)
        x0 = Float64(x0r)
        if x0 <= x[1]
            j = 1
        elseif x0 >= x[end]
            j = n - 1
        else
            j = searchsortedlast(x, x0)
            j = clamp(j, 1, n - 1)
        end
        t = (x0 - x[j]) / (x[j+1] - x[j])
        out[i] = (1.0 - t)*y[j] + t*y[j+1]
    end
    return out
end

function make_interpolated_background(Egrid::AbstractVector{<:Real},
                                      Eraw::AbstractVector{<:Real},
                                      Iraw::AbstractVector{<:Real};
                                      smooth_sigma_meV::Real=0.0,
                                      interpolation_kind::Symbol=:pchip)
    Es, Is = sort_xy(Eraw, Iraw)
    Is_smooth = gaussian_smooth_xy(Es, Is, smooth_sigma_meV)

    bg = if interpolation_kind == :linear
        linear_interpolate(Es, Is_smooth, Egrid)
    elseif interpolation_kind == :pchip
        # Shape-preserving cubic Hermite interpolation.  This is smooth, but unlike a
        # natural cubic spline it strongly suppresses overshoot and wild curvature.
        pchip_interpolator(Es, Is_smooth)(Egrid)
    elseif interpolation_kind == :spline
        # Kept only for comparison/debugging; not recommended for these backgrounds.
        natural_cubic_spline(Es, Is_smooth)(Egrid)
    else
        error("Unknown interpolation_kind=$(interpolation_kind). Use :pchip, :linear, or :spline.")
    end

    return bg, Es, Is_smooth
end


# -------------------------------------------------------------------------
# Quasi-physics continuum fits for the 0 T baseline under the structured peak
# -------------------------------------------------------------------------

function ridge_linear_least_squares(X::AbstractMatrix{<:Real},
                                    y::AbstractVector{<:Real},
                                    err::Union{Nothing,AbstractVector{<:Real}}=nothing;
                                    ridge_lambda::Real=0.0,
                                    unpenalized_columns::Vector{Int}=Int[])
    Xf = Float64.(X)
    yf = Float64.(y)
    if err !== nothing
        σ = max.(Float64.(err), eps(Float64))
        sw = 1.0 ./ σ
        Xf = Xf .* sw
        yf = yf .* sw
    end

    A = transpose(Xf) * Xf
    b = transpose(Xf) * yf

    if ridge_lambda > 0
        penalty = Matrix{Float64}(I, size(A, 1), size(A, 2))
        for j in unpenalized_columns
            penalty[j, j] = 0.0
        end
        A .+= Float64(ridge_lambda) .* penalty
    end

    return A \ b
end

function continuum_design_matrix(E::AbstractVector{<:Real};
                                 model::Symbol=:power_tail,
                                 power::Real=3.0,
                                 energy_offset::Real=0.15,
                                 include_linear_tilt::Bool=false)
    Ef = Float64.(E)
    cols = Vector{Vector{Float64}}()
    push!(cols, ones(length(Ef)))

    if model == :power_tail
        push!(cols, (Ef .+ Float64(energy_offset)) .^ (-Float64(power)))
    elseif model == :exp_tail
        # Here `power` is used as tau in meV for the exponential tail.
        push!(cols, exp.(-Ef ./ Float64(power)))
    elseif model == :inverse_series
        push!(cols, (Ef .+ Float64(energy_offset)) .^ (-1.0))
        push!(cols, (Ef .+ Float64(energy_offset)) .^ (-2.0))
    elseif model == :line
        push!(cols, Ef .- mean(Ef))
    else
        error("Unknown continuum model $(model). Use :power_tail, :exp_tail, :inverse_series, or :line.")
    end

    if include_linear_tilt && model != :line
        push!(cols, Ef .- mean(Ef))
    end

    X = zeros(Float64, length(Ef), length(cols))
    for j in eachindex(cols)
        X[:, j] .= cols[j]
    end
    return X
end

function eval_continuum_model(E::AbstractVector{<:Real}, coeff::AbstractVector{<:Real};
                              model::Symbol=:power_tail,
                              power::Real=3.0,
                              energy_offset::Real=0.15,
                              include_linear_tilt::Bool=false,
                              center_energy::Real=0.0)
    Ef = Float64.(E)
    y = coeff[1] .* ones(length(Ef))
    j = 2

    if model == :power_tail
        y .+= coeff[j] .* (Ef .+ Float64(energy_offset)) .^ (-Float64(power))
        j += 1
    elseif model == :exp_tail
        y .+= coeff[j] .* exp.(-Ef ./ Float64(power))
        j += 1
    elseif model == :inverse_series
        y .+= coeff[j] .* (Ef .+ Float64(energy_offset)) .^ (-1.0)
        j += 1
        y .+= coeff[j] .* (Ef .+ Float64(energy_offset)) .^ (-2.0)
        j += 1
    elseif model == :line
        y .+= coeff[j] .* (Ef .- Float64(center_energy))
        j += 1
    else
        error("Unknown continuum model $(model)")
    end

    if include_linear_tilt && model != :line
        y .+= coeff[j] .* (Ef .- Float64(center_energy))
    end
    return y
end

function fit_zeroT_continuum_baseline(Efit::AbstractVector{<:Real},
                                      Ifit::AbstractVector{<:Real};
                                      errfit::Union{Nothing,AbstractVector{<:Real}}=nothing,
                                      model::Symbol=:power_tail,
                                      power_grid=collect(0.5:0.05:8.0),
                                      exp_tau_grid=collect(0.25:0.025:3.0),
                                      energy_offset::Real=0.15,
                                      ridge_lambda::Real=0.0,
                                      include_linear_tilt::Bool=false,
                                      positive_tail::Bool=true)
    Ef = Float64.(Efit)
    If = Float64.(Ifit)
    center_energy = mean(Ef)

    scan_grid = if model == :power_tail
        collect(power_grid)
    elseif model == :exp_tail
        collect(exp_tau_grid)
    else
        [NaN]
    end

    best_score = Inf
    best_coeff = Float64[]
    best_param = NaN

    for param in scan_grid
        X = continuum_design_matrix(Ef;
            model=model,
            power=isnan(param) ? 3.0 : param,
            energy_offset=energy_offset,
            include_linear_tilt=include_linear_tilt,
        )
        coeff = ridge_linear_least_squares(X, If, errfit;
            ridge_lambda=ridge_lambda,
            unpenalized_columns=[1],
        )

        # For the tail models, coeff[2] is the decaying tail amplitude.  Keeping it
        # positive prevents the baseline from turning into an unphysical rising tail.
        if positive_tail && (model == :power_tail || model == :exp_tail) && coeff[2] < 0
            continue
        end

        pred = X * coeff
        resid = If .- pred
        if errfit !== nothing
            σ = max.(Float64.(errfit), eps(Float64))
            score = mean((resid ./ σ).^2)
        else
            score = mean(resid.^2)
        end

        if score < best_score
            best_score = score
            best_coeff = Float64.(coeff)
            best_param = isnan(param) ? 3.0 : Float64(param)
        end
    end

    if isempty(best_coeff)
        # Fallback: allow the sign if every positive-tail fit was rejected.
        return fit_zeroT_continuum_baseline(Efit, Ifit;
            errfit=errfit,
            model=model,
            power_grid=power_grid,
            exp_tau_grid=exp_tau_grid,
            energy_offset=energy_offset,
            ridge_lambda=ridge_lambda,
            include_linear_tilt=include_linear_tilt,
            positive_tail=false,
        )
    end

    baseline(Enew) = eval_continuum_model(Enew, best_coeff;
        model=model,
        power=best_param,
        energy_offset=energy_offset,
        include_linear_tilt=include_linear_tilt,
        center_energy=center_energy,
    )

    return baseline, best_coeff, best_param, best_score
end

function structured_residual_points(byfield::Dict{Float64,Scan1D},
                                    residual_window::Tuple{Float64,Float64};
                                    fit_window::Tuple{Float64,Float64}=STRUCTURED_FIT_WINDOW,
                                    model::Symbol=:power_tail,
                                    power_grid=collect(0.5:0.05:8.0),
                                    exp_tau_grid=collect(0.25:0.025:3.0),
                                    energy_offset::Real=0.15,
                                    ridge_lambda::Real=0.0,
                                    include_linear_tilt::Bool=false,
                                    positive_tail::Bool=true,
                                    clip_negative_residuals::Bool=false)
    s0 = byfield[0.0]
    E = s0.energy
    I0 = s0.intensity
    fit_lo, fit_hi = fit_window
    res_lo, res_hi = residual_window

    peak_mask = (E .>= res_lo) .& (E .<= res_hi)
    fit_mask = (E .>= fit_lo) .& (E .<= fit_hi) .& .!peak_mask

    any(fit_mask) || error("No points found for 0 T continuum fit outside residual window")
    any(peak_mask) || error("No points found inside residual window $(residual_window)")

    Efit = E[fit_mask]
    Ifit = I0[fit_mask]
    errfit = s0.error[fit_mask]

    baseline_fun, coeff, param, score = fit_zeroT_continuum_baseline(Efit, Ifit;
        errfit=errfit,
        model=model,
        power_grid=power_grid,
        exp_tau_grid=exp_tau_grid,
        energy_offset=energy_offset,
        ridge_lambda=ridge_lambda,
        include_linear_tilt=include_linear_tilt,
        positive_tail=positive_tail,
    )

    Eres = E[peak_mask]
    baseline = baseline_fun(Eres)
    residual = I0[peak_mask] .- baseline
    if clip_negative_residuals
        residual = max.(residual, 0.0)
    end

    @printf("%s 0T continuum baseline: model=%s, param=%.4g, score=%.4g, coeff=%s\n",
            s0.qtag, String(model), param, score, repr(coeff))

    return Eres, residual, Efit, Ifit, baseline
end


function make_spline_background_model(qtag::String,
                                      byfield::Dict{Float64,Scan1D};
                                      fields=FIELD_ORDER,
                                      low_window=MIN_BG_LOW_WINDOW,
                                      high_threshold=MIN_BG_HIGH_THRESHOLD,
                                      residual_windows=STRUCTURED_RESIDUAL_WINDOWS,
                                      final_smooth_sigma_meV::Real=0.0,
                                      final_interp_kind::Symbol=:pchip,
                                      fit_window::Tuple{Float64,Float64}=STRUCTURED_FIT_WINDOW,
                                      zeroT_baseline_model::Symbol=:power_tail,
                                      power_grid=collect(0.5:0.05:8.0),
                                      exp_tau_grid=collect(0.25:0.025:3.0),
                                      energy_offset::Real=0.15,
                                      ridge_lambda::Real=0.0,
                                      include_linear_tilt::Bool=false,
                                      positive_tail::Bool=true,
                                      clip_negative_residuals::Bool=false)
    Egrid = common_energy_grid(byfield; fields=fields)
    E_lowhigh, I_lowhigh = min_over_fields_background_raw(byfield;
        fields=fields,
        low_window=low_window,
        high_threshold=high_threshold,
    )

    Eraw = copy(E_lowhigh)
    Iraw = copy(I_lowhigh)

    residual_window = haskey(residual_windows, qtag) ? residual_windows[qtag] : nothing
    Eres = Float64[]
    residual = Float64[]
    residual_abs_bg = Float64[]
    Efit = Float64[]
    Ifit = Float64[]
    baseline_reswin = Float64[]

    if residual_window !== nothing
        Eres, residual, Efit, Ifit, baseline_reswin = structured_residual_points(byfield, residual_window;
            fit_window=fit_window,
            model=zeroT_baseline_model,
            power_grid=power_grid,
            exp_tau_grid=exp_tau_grid,
            energy_offset=energy_offset,
            ridge_lambda=ridge_lambda,
            include_linear_tilt=include_linear_tilt,
            positive_tail=positive_tail,
            clip_negative_residuals=clip_negative_residuals,
        )

        # Convert the residual peak into absolute background points by placing it on top
        # of the smooth low/high minimum-field background estimate.
        bg_without_residual, _, _ = make_interpolated_background(Egrid, E_lowhigh, I_lowhigh;
            smooth_sigma_meV=final_smooth_sigma_meV,
            interpolation_kind=final_interp_kind,
        )
        residual_base = linear_interpolate(Float64.(Egrid), Float64.(bg_without_residual), Eres)
        residual_abs_bg = residual_base .+ residual

        append!(Eraw, Eres)
        append!(Iraw, residual_abs_bg)
    end

    bg, Esorted, Ismoothed = make_interpolated_background(Egrid, Eraw, Iraw;
        smooth_sigma_meV=final_smooth_sigma_meV,
        interpolation_kind=final_interp_kind,
    )

    return BackgroundModel(
        qtag,
        Esorted,
        Float64.(Iraw[sortperm(Eraw)]),
        Ismoothed,
        Float64.(Egrid),
        Float64.(bg),
        Float64.(E_lowhigh),
        Float64.(I_lowhigh),
        residual_window,
        Float64.(Eres),
        Float64.(residual),
        Float64.(residual_abs_bg),
        Float64.(Efit),
        Float64.(Ifit),
        Float64.(baseline_reswin),
    )
end

function subtract_background(scan::Scan1D, bg::BackgroundModel)
    length(scan.energy) == length(bg.energy_grid) || error("Background length mismatch for $(scan.filename)")
    maximum(abs.(scan.energy .- bg.energy_grid)) <= 1e-10 || error("Energy grid mismatch for $(scan.filename)")

    return Scan1D(
        scan.qtag,
        scan.field_T,
        scan.temperature_K,
        scan.Ei_meV,
        scan.filename,
        scan.meta,
        scan.intensity .- bg.background,
        copy(scan.error),
        copy(scan.energy),
        copy(scan.K),
        copy(scan.L),
        copy(scan.H),
    )
end

function make_spline_background_subtracted_scans(scans::Dict{String,Dict{Float64,Scan1D}};
                                                 kwargs...)
    corrected = Dict{String,Dict{Float64,Scan1D}}()
    background_models = Dict{String,BackgroundModel}()

    for qtag in Q_ORDER
        haskey(scans, qtag) || continue
        byfield = scans[qtag]
        bg = make_spline_background_model(qtag, byfield; kwargs...)
        background_models[qtag] = bg

        corrected[qtag] = Dict{Float64,Scan1D}()
        for (B, s) in byfield
            corrected[qtag][B] = subtract_background(s, bg)
        end
    end

    return corrected, background_models
end

function plot_raw_and_spline_bgsubtracted(scans::Dict{String,Dict{Float64,Scan1D}};
                                          outdir::AbstractString=OUTDIR,
                                          save_png::Bool=true,
                                          raw_ylim::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                          bgsub_ylim::Tuple{Float64,Float64}=(-0.001, 0.004),
                                          final_smooth_sigma_meV::Real=0.0,
                                          final_interp_kind::Symbol=:pchip,
                                          zeroT_baseline_model::Symbol=:power_tail,
                                          power_grid=collect(0.5:0.05:8.0),
                                          exp_tau_grid=collect(0.25:0.025:3.0),
                                          energy_offset::Real=0.15,
                                          ridge_lambda::Real=0.0,
                                          include_linear_tilt::Bool=false,
                                          positive_tail::Bool=true,
                                          clip_negative_residuals::Bool=false)
    mkpath(outdir)

    corrected, background_models = make_spline_background_subtracted_scans(scans;
        final_smooth_sigma_meV=final_smooth_sigma_meV,
        final_interp_kind=final_interp_kind,
        zeroT_baseline_model=zeroT_baseline_model,
        power_grid=power_grid,
        exp_tau_grid=exp_tau_grid,
        energy_offset=energy_offset,
        ridge_lambda=ridge_lambda,
        include_linear_tilt=include_linear_tilt,
        positive_tail=positive_tail,
        clip_negative_residuals=clip_negative_residuals,
    )

    fig = Figure(size=(1450, 850))

    for (icol, qtag) in enumerate(Q_ORDER)
        haskey(scans, qtag) || continue
        meta = QMETA[qtag]
        ylim_top = raw_ylim === nothing ? meta.plot_ylim : raw_ylim

        ax_raw = Axis(fig[1, icol];
            title=meta.label,
            xlabel=icol == 2 ? "ΔE (meV)" : "",
            ylabel=icol == 1 ? "Raw intensity" : "",
            limits=(meta.de_range[1], meta.de_range[2], ylim_top[1], ylim_top[2]),
        )

        for B in sort(collect(keys(scans[qtag])))
            s = scans[qtag][B]
            lines!(ax_raw, s.energy, s.intensity; label=@sprintf("%.0f T", B), linewidth=2)
            errorbars!(ax_raw, s.energy, s.intensity, s.error; whiskerwidth=4)
        end
        axislegend(ax_raw; position=:rt, framevisible=false)

        ax_sub = Axis(fig[2, icol];
            xlabel="ΔE (meV)",
            ylabel=icol == 1 ? "Intensity - background" : "",
            limits=(meta.de_range[1], meta.de_range[2], bgsub_ylim[1], bgsub_ylim[2]),
        )

        for B in sort(collect(keys(corrected[qtag])))
            s = corrected[qtag][B]
            lines!(ax_sub, s.energy, s.intensity; label=@sprintf("%.0f T", B), linewidth=2)
            errorbars!(ax_sub, s.energy, s.intensity, s.error; whiskerwidth=4)
        end

        bg = background_models[qtag]
        lines!(ax_sub, bg.energy_grid, zeros(length(bg.energy_grid)); linestyle=:dash, linewidth=1)
    end

    Label(fig[0, :], "YZGO 1D energy scans: raw and model-background-subtracted", fontsize=20)

    if save_png
        save(joinpath(outdir, "YZGO_1d_Escan_raw_and_tail_bgsubtracted.png"), fig)
    end

    return fig, corrected, background_models
end


# New clearer name. The old function name is retained for backward compatibility.
plot_raw_and_tail_bgsubtracted(scans::Dict{String,Dict{Float64,Scan1D}}; kwargs...) =
    plot_raw_and_spline_bgsubtracted(scans; kwargs...)

function plot_tail_background_models(background_models::Dict{String,BackgroundModel}; kwargs...)
    return plot_spline_background_models(background_models; kwargs...)
end

function plot_zeroT_tail_residual_diagnostics(scans::Dict{String,Dict{Float64,Scan1D}},
                                              background_models::Dict{String,BackgroundModel};
                                              kwargs...)
    return plot_zeroT_spline_residual_diagnostics(scans, background_models; kwargs...)
end

function plot_spline_background_models(background_models::Dict{String,BackgroundModel};
                                       outdir::AbstractString=OUTDIR,
                                       save_png::Bool=true)
    mkpath(outdir)
    fig = Figure(size=(1450, 430))

    for (icol, qtag) in enumerate(Q_ORDER)
        haskey(background_models, qtag) || continue
        bg = background_models[qtag]
        meta = QMETA[qtag]

        ax = Axis(fig[1, icol];
            title=meta.label,
            xlabel="ΔE (meV)",
            ylabel=icol == 1 ? "Background estimate" : "",
        )
        xlims!(ax, meta.de_range[1], meta.de_range[2])

        scatter!(ax, bg.lowhigh_energy, bg.lowhigh_intensity; label="min over fields anchors", markersize=7)
        if bg.residual_window !== nothing
            scatter!(ax, bg.residual_energy, bg.residual_absolute_background;
                label="structured residual anchors", markersize=9)
        end
        lines!(ax, bg.energy_grid, bg.background; label="final shape-preserving background", linewidth=3)
        axislegend(ax; position=:rt, framevisible=false)
    end

    if save_png
        save(joinpath(outdir, "YZGO_1d_Escan_tail_background_models.png"), fig)
    end
    return fig
end

function plot_zeroT_spline_residual_diagnostics(scans::Dict{String,Dict{Float64,Scan1D}},
                                                background_models::Dict{String,BackgroundModel};
                                                outdir::AbstractString=OUTDIR,
                                                save_png::Bool=true)
    mkpath(outdir)
    qtags = [q for q in Q_ORDER if haskey(background_models, q) && background_models[q].residual_window !== nothing]
    fig = Figure(size=(1000, 430))

    for (icol, qtag) in enumerate(qtags)
        bg = background_models[qtag]
        s0 = scans[qtag][0.0]
        meta = QMETA[qtag]
        ax = Axis(fig[1, icol];
            title=meta.label,
            xlabel="ΔE (meV)",
            ylabel=icol == 1 ? "0 T intensity / continuum fit" : "",
        )
        xlims!(ax, 0.9, 3.1)

        lines!(ax, s0.energy, s0.intensity; label="0 T data", linewidth=2)
        scatter!(ax, bg.zeroT_fit_energy, bg.zeroT_fit_intensity; label="fit points", markersize=6)
        lines!(ax, bg.residual_energy, bg.zeroT_baseline_in_residual_window;
            label="continuum baseline across excluded window", linewidth=3)
        axislegend(ax; position=:rt, framevisible=false)
    end

    if save_png
        save(joinpath(outdir, "YZGO_1d_Escan_zeroT_tail_residual_diagnostics.png"), fig)
    end
    return fig
end

# Backward-compatible simple constant-background tools from the first version.
function constant_background(scan::Scan1D; windows::Vector{Tuple{Float64,Float64}})
    mask = energy_window_mask(scan.energy, windows)
    any(mask) || error("No points found in windows $(windows) for $(scan.filename)")

    σ = max.(scan.error[mask], eps(Float64))
    w = 1.0 ./ σ.^2
    bg = sum(w .* scan.intensity[mask]) / sum(w)
    bgerr = sqrt(1.0 / sum(w))
    return bg, bgerr, mask
end

function subtract_constant_background(scan::Scan1D; windows::Vector{Tuple{Float64,Float64}})
    bg, bgerr, _ = constant_background(scan; windows=windows)
    return Scan1D(
        scan.qtag,
        scan.field_T,
        scan.temperature_K,
        scan.Ei_meV,
        scan.filename,
        scan.meta,
        scan.intensity .- bg,
        sqrt.(scan.error.^2 .+ bgerr^2),
        copy(scan.energy),
        copy(scan.K),
        copy(scan.L),
        copy(scan.H),
    )
end

function make_bg_subtracted_scans(scans::Dict{String,Dict{Float64,Scan1D}};
                                  windows_by_q::Dict{String,Vector{Tuple{Float64,Float64}}})
    corrected = Dict{String,Dict{Float64,Scan1D}}()
    bg_values = Dict{Tuple{String,Float64},Tuple{Float64,Float64}}()

    for (qtag, byfield) in scans
        windows = windows_by_q[qtag]
        corrected[qtag] = Dict{Float64,Scan1D}()
        for (B, s) in byfield
            bg, bgerr, _ = constant_background(s; windows=windows)
            bg_values[(qtag, B)] = (bg, bgerr)
            corrected[qtag][B] = subtract_constant_background(s; windows=windows)
        end
    end
    return corrected, bg_values
end

const BG_WINDOWS_EXAMPLE = Dict(
    "0_1_0" => [(2.8, 3.25)],
    "0p33_0p33_0" => [(3.4, 3.95)],
    "0p5_0_0" => [(3.4, 3.95)],
)

function plot_bg_subtracted_example(scans::Dict{String,Dict{Float64,Scan1D}};
                                    windows_by_q=BG_WINDOWS_EXAMPLE,
                                    outdir::AbstractString=OUTDIR)
    corrected, bg_values = make_bg_subtracted_scans(scans; windows_by_q=windows_by_q)
    fig = plot_initial_scans(corrected; outdir=outdir, save_png=false)
    save(joinpath(outdir, "YZGO_1d_Escan_bg_subtracted_example.png"), fig)
    return fig, corrected, bg_values
end


# =============================================================================
# Section copied from YZGO_analytic_1d2d_glmakie_stratified_1d2d.jl, with auto-run block removed
# =============================================================================

# YZGO_analytic_1d2d_glmakie.jl
#
# Shared analytical field-polarized triangular-lattice J1-J2 model for YZGO,
# with optional disorder, optional instrumental energy resolution, optional
# momentum/mosaic broadening, and GLMakie plotting for both 1D cuts and a 2D
# E-vs-momentum map along a chosen reciprocal-space leg.
#
# Requires: SpecialFunctions, GLMakie
#
# Typical usage:
#     include("YZGO_analytic_1d2d_glmakie.jl")
#     run_demo_1d()
#     run_demo_2d()
#
# At the bottom of this file, the script entry point is commented so you can
# choose whether to run the 1D or 2D demo when executing as a script.

using Random
using LinearAlgebra
using Printf
using SpecialFunctions
using GLMakie

# -----------------------------
# Constants
# -----------------------------

const MU_B_MEV_PER_T = 5.7883818060e-2
const FWHM_TO_SIGMA = 1.0 / (2.0 * sqrt(2.0 * log(2.0)))
const SIGMA_TO_FWHM = 1.0 / FWHM_TO_SIGMA

# -----------------------------
# Parameter containers
# -----------------------------

Base.@kwdef struct ModelParams
    B_T::Float64
    gzz::Float64
    J1_meV::Float64
    J2_meV::Float64
    S::Float64 = 0.5
end

Base.@kwdef struct DisorderParams
    enabled::Bool = true
    sigma_gzz::Float64 = 1.0 / 6.0
    sigma_J1::Float64 = 2.0 / 3.0
    sigma_J2::Float64 = 2.0 / 3.0
    J_units::Symbol = :fractional   # :fractional or :absolute_meV
    correlate_J1_J2::Bool = false
end

Base.@kwdef struct ResolutionParams
    Ei_meV::Float64 = 4.65
    energy_enabled::Bool = true

    # Momentum broadening model.
    #   :rlu_diagonal   -> Gaussian in H,K,L directly.
    #   :cartesian_Ainv -> isotropic Gaussian in |Q| space.
    momentum_enabled::Bool = true
    momentum_mode::Symbol = :rlu_diagonal

    # Mosaic / resolution widths in reciprocal lattice units.
    sigma_H_rlu::Float64 = 0.037
    sigma_K_rlu::Float64 = 0.037
    sigma_L_rlu::Float64 = 0.252

    # Retained for comparisons with the old isotropic-|Q| approximation.
    sigma_Q_Ainv::Float64 = 0.017

    # Expand the TRUE-Q sampling region by this many sigmas outside the
    # nominal experimental histogram / integration bounds to capture bleed-in.
    sample_margin_nsigma::Float64 = 3.0
end

Base.@kwdef struct LatticeParams
    a_A::Float64 = 3.376
    c_A::Float64 = 21.96
end

Base.@kwdef struct CutSpec1D
    name::String
    H_range::Tuple{Float64, Float64}
    K_range::Tuple{Float64, Float64}
    L_range::Tuple{Float64, Float64}
    E_range::Tuple{Float64, Float64}
    dE_meV::Float64 = 0.05
end

Base.@kwdef struct PathCutSpec2D
    name::String
    uvec::NTuple{3, Float64}
    vvec::NTuple{3, Float64}
    wvec::NTuple{3, Float64}
    u_range::Tuple{Float64, Float64}
    v_range::Tuple{Float64, Float64}
    w_range::Tuple{Float64, Float64}
    du::Float64 = 0.01
    E_range::Tuple{Float64, Float64} = (0.0, 4.0)
    dE_meV::Float64 = 0.05
end

# -----------------------------
# Default experimental cut specs
# -----------------------------

function default_cuts_1d()
    return [
        CutSpec1D(
            name = "Gamma_cut_center_0_1_0",
            H_range = (-0.1, 0.1),
            K_range = (0.9, 1.1),
            L_range = (-0.3, 0.3),
            E_range = (0.0, 3.3),
            dE_meV = 0.05,
        ),
        # Experimental label M, but in the simple triangular-lattice
        # convention used here, this is the analytic K-point region.
        CutSpec1D(
            name = "M_label_cut_center_1over3_1over3_0",
            H_range = (0.23, 0.43),
            K_range = (0.23, 0.43),
            L_range = (-0.3, 0.3),
            E_range = (0.0, 4.0),
            dE_meV = 0.05,
        ),
        # Experimental label K, but in the simple triangular-lattice
        # convention used here, this is the analytic M-point region.
        CutSpec1D(
            name = "K_label_cut_center_1over2_0_0",
            H_range = (0.4, 0.6),
            K_range = (-0.1, 0.1),
            L_range = (-0.3, 0.3),
            E_range = (0.0, 4.0),
            dE_meV = 0.05,
        ),
    ]
end

function default_leg_cut_2d()
    return PathCutSpec2D(
        name = "K1_to_Gamma1_leg",
        uvec = (1.0, -0.5, 0.0),
        vvec = (0.0, 1.0, 0.0),
        wvec = (0.0, 0.0, 1.0),
        u_range = (-1.0 / 3.0, 1.0),
        v_range = (0.45, 0.55),
        w_range = (-0.25, 0.25),
        du = 0.01,
        E_range = (0.0, 4.0),
        dE_meV = 0.05,
    )
end

# -----------------------------
# Reciprocal-space conversions
# -----------------------------

function rlu_basis_matrix(lattice::LatticeParams)
    astar = 4.0 * pi / (sqrt(3.0) * lattice.a_A)
    cstar = 2.0 * pi / lattice.c_A

    b1 = [ astar, 0.0, 0.0 ]
    b2 = [ -0.5 * astar, 0.5 * sqrt(3.0) * astar, 0.0 ]
    b3 = [ 0.0, 0.0, cstar ]

    return hcat(b1, b2, b3)
end

function rlu_to_cart(H::Real, K::Real, L::Real, lattice::LatticeParams)
    Bmat = rlu_basis_matrix(lattice)
    return Bmat * [Float64(H), Float64(K), Float64(L)]
end

function cart_to_rlu(qcart::AbstractVector{<:Real}, lattice::LatticeParams)
    Bmat = rlu_basis_matrix(lattice)
    return Bmat \ collect(Float64, qcart)
end

# -----------------------------
# 2D path coordinate transforms
# -----------------------------

function projection_matrix(spec::PathCutSpec2D)
    u = collect(spec.uvec)
    v = collect(spec.vvec)
    w = collect(spec.wvec)
    return hcat(u, v, w)
end

function uvw_to_hkl(u::Real, v::Real, w::Real, spec::PathCutSpec2D)
    P = projection_matrix(spec)
    hkl = P * [Float64(u), Float64(v), Float64(w)]
    return (hkl[1], hkl[2], hkl[3])
end

function hkl_to_uvw(H::Real, K::Real, L::Real, spec::PathCutSpec2D)
    P = projection_matrix(spec)
    uvw = P \ [Float64(H), Float64(K), Float64(L)]
    return (uvw[1], uvw[2], uvw[3])
end

# -----------------------------
# Analytical dispersion
# -----------------------------

function Delta1_triangular(H::Real, K::Real)
    return 6.0 - 2.0 * (
        cos(2.0 * pi * H) +
        cos(2.0 * pi * K) +
        cos(2.0 * pi * (H + K))
    )
end

function Delta2_triangular(H::Real, K::Real)
    return 6.0 - 2.0 * (
        cos(2.0 * pi * (H - K)) +
        cos(2.0 * pi * (2.0 * H + K)) +
        cos(2.0 * pi * (H + 2.0 * K))
    )
end

function dispersion_meV(H::Real, K::Real, model::ModelParams;
                        gzz::Real = model.gzz,
                        J1_meV::Real = model.J1_meV,
                        J2_meV::Real = model.J2_meV)
    return gzz * MU_B_MEV_PER_T * model.B_T -
           model.S * (J1_meV * Delta1_triangular(H, K) +
                      J2_meV * Delta2_triangular(H, K))
end

function fit_model_from_GKM(EGamma_meV::Real, E_Kanalytic_meV::Real, E_Manalytic_meV::Real, B_T::Real;
                            S::Real = 0.5)
    EGamma = Float64(EGamma_meV)
    EK = Float64(E_Kanalytic_meV)
    EM = Float64(E_Manalytic_meV)
    B = Float64(B_T)
    Sval = Float64(S)

    gzz = EGamma / (MU_B_MEV_PER_T * B)
    J1 = (EGamma - EK) / (9.0 * Sval)
    J2 = (EGamma - EM) / (8.0 * Sval) - J1

    return ModelParams(B_T = B, gzz = gzz, J1_meV = J1, J2_meV = J2, S = Sval)
end

# -----------------------------
# Energy resolution
# -----------------------------

function energy_resolution_fwhm_meV(Etransfer_meV::Real, Ei_meV::Real = 4.65)
    E = Float64(Etransfer_meV)
    Ei = Float64(Ei_meV)

    if !isfinite(E) || E < 0.0 || E >= Ei
        return NaN
    end

    ef = Ei - E
    if ef <= 0.0
        return NaN
    end

    return 2.4994e-4 * sqrt(
        ef^3 * (
            (146.25936 * (0.052 + 0.123 * (Ei / ef)^1.5))^2 +
            (57.27700  * (1.052 + 0.123 * (Ei / ef)^1.5))^2
        )
    )
end

energy_resolution_sigma_meV(Etransfer_meV::Real, Ei_meV::Real = 4.65) =
    FWHM_TO_SIGMA * energy_resolution_fwhm_meV(Etransfer_meV, Ei_meV)

# -----------------------------
# Yb3+ magnetic form factor
# -----------------------------

function radial_integral_j0(s::Real, A::Real, a::Real, B::Real, b::Real, C::Real, c::Real, D::Real)
    x = Float64(s)^2
    return A * exp(-a * x) + B * exp(-b * x) + C * exp(-c * x) + D
end

function radial_integral_jL_nonzero(s::Real, A::Real, a::Real, B::Real, b::Real, C::Real, c::Real, D::Real)
    x = Float64(s)^2
    return (A * exp(-a * x) + B * exp(-b * x) + C * exp(-c * x) + D) * x
end

function yb3_j0(Q_Ainv::Real)
    s = Float64(Q_Ainv) / (4.0 * pi)
    return radial_integral_j0(s, 0.0416, 16.0949, 0.2849, 7.8341, 0.6961, 2.6725, -0.0229)
end

function yb3_j2(Q_Ainv::Real)
    s = Float64(Q_Ainv) / (4.0 * pi)
    return radial_integral_jL_nonzero(s, 0.1570, 18.5553, 0.8484, 6.5403, 0.8880, 2.0367, 0.0318)
end

function yb3_form_factor(Q_Ainv::Real; include_j2::Bool = true, c2::Real = 0.75)
    if include_j2
        return yb3_j0(Q_Ainv) + Float64(c2) * yb3_j2(Q_Ainv)
    else
        return yb3_j0(Q_Ainv)
    end
end

# -----------------------------
# Cross-section weights
# -----------------------------

function polarization_factor(qcart::AbstractVector{<:Real}; mode::Symbol = :transverse_c)
    Q = norm(qcart)
    if Q < 1e-12 || mode == :none
        return 1.0
    end

    qzhat = Float64(qcart[3]) / Q

    if mode == :transverse_c
        return 0.5 * (1.0 + qzhat^2)
    elseif mode == :static_c
        return max(0.0, 1.0 - qzhat^2)
    else
        error("Unknown polarization mode: $mode. Use :transverse_c, :static_c, or :none.")
    end
end

function kf_over_ki(Etransfer_meV::Real, Ei_meV::Real)
    E = Float64(Etransfer_meV)
    Ei = Float64(Ei_meV)
    Ef = Ei - E
    if Ef <= 0.0
        return 0.0
    end
    return sqrt(Ef / Ei)
end

# -----------------------------
# Histogram utilities
# -----------------------------

bin_centers(edges::AbstractVector{<:Real}) = 0.5 .* (edges[1:end-1] .+ edges[2:end])
rand_uniform(rng::AbstractRNG, range::Tuple{Float64, Float64}) = range[1] + rand(rng) * (range[2] - range[1])
in_range(x::Real, range::Tuple{Float64, Float64}) = range[1] <= Float64(x) <= range[2]

function gaussian_bin_probability(lo::Real, hi::Real, mu::Real, sigma::Real)
    if sigma <= 0.0 || !isfinite(sigma)
        return (lo <= mu < hi) ? 1.0 : 0.0
    end
    zlo = (Float64(lo) - Float64(mu)) / (sqrt(2.0) * Float64(sigma))
    zhi = (Float64(hi) - Float64(mu)) / (sqrt(2.0) * Float64(sigma))
    return 0.5 * (erf(zhi) - erf(zlo))
end

function deposit_delta_1d!(hist::Vector{Float64}, edges::Vector{Float64}, E::Real, weight::Real)
    idx = searchsortedlast(edges, Float64(E))
    if idx >= 1 && idx <= length(hist)
        hist[idx] += Float64(weight)
    end
    return nothing
end

function deposit_gaussian_1d!(hist::Vector{Float64}, edges::Vector{Float64}, mu::Real, sigma::Real, weight::Real;
                              nsigma_window::Real = 6.0)
    sigma = Float64(sigma)
    mu = Float64(mu)
    weight = Float64(weight)

    if sigma <= 0.0 || !isfinite(sigma)
        deposit_delta_1d!(hist, edges, mu, weight)
        return nothing
    end

    loE = mu - nsigma_window * sigma
    hiE = mu + nsigma_window * sigma

    i0 = max(1, searchsortedlast(edges, loE))
    i1 = min(length(hist), searchsortedfirst(edges, hiE))

    for i in i0:i1
        pbin = gaussian_bin_probability(edges[i], edges[i + 1], mu, sigma)
        hist[i] += weight * pbin
    end
    return nothing
end

function deposit_2d_energyline!(hist::Matrix{Float64}, uedges::Vector{Float64}, eedges::Vector{Float64},
                                umeas::Real, emu::Real, esigma::Real, weight::Real;
                                nsigma_window::Real = 6.0)
    iu = searchsortedlast(uedges, Float64(umeas))
    if !(1 <= iu <= size(hist, 1))
        return nothing
    end

    sigma = Float64(esigma)
    mu = Float64(emu)
    weight = Float64(weight)

    if sigma <= 0.0 || !isfinite(sigma)
        ie = searchsortedlast(eedges, mu)
        if 1 <= ie <= size(hist, 2)
            hist[iu, ie] += weight
        end
        return nothing
    end

    loE = mu - nsigma_window * sigma
    hiE = mu + nsigma_window * sigma
    j0 = max(1, searchsortedlast(eedges, loE))
    j1 = min(size(hist, 2), searchsortedfirst(eedges, hiE))

    for j in j0:j1
        pbin = gaussian_bin_probability(eedges[j], eedges[j + 1], mu, sigma)
        hist[iu, j] += weight * pbin
    end
    return nothing
end

# -----------------------------
# Momentum resolution / mosaic helpers
# -----------------------------

function momentum_sigmas_hkl(resolution::ResolutionParams, lattice::LatticeParams)
    if !resolution.momentum_enabled
        return (0.0, 0.0, 0.0)
    end

    if resolution.momentum_mode == :rlu_diagonal
        return (
            resolution.sigma_H_rlu,
            resolution.sigma_K_rlu,
            resolution.sigma_L_rlu,
        )
    elseif resolution.momentum_mode == :cartesian_Ainv
        Bmat = rlu_basis_matrix(lattice)
        return (
            resolution.sigma_Q_Ainv / norm(Bmat[:, 1]),
            resolution.sigma_Q_Ainv / norm(Bmat[:, 2]),
            resolution.sigma_Q_Ainv / norm(Bmat[:, 3]),
        )
    else
        error("Unknown momentum_mode = $(resolution.momentum_mode). Use :rlu_diagonal or :cartesian_Ainv.")
    end
end

function smear_hkl_to_measured(Htrue::Real, Ktrue::Real, Ltrue::Real,
                               resolution::ResolutionParams,
                               lattice::LatticeParams,
                               rng::AbstractRNG)
    if !resolution.momentum_enabled
        return (Float64(Htrue), Float64(Ktrue), Float64(Ltrue))
    end

    if resolution.momentum_mode == :rlu_diagonal
        return (
            Float64(Htrue) + resolution.sigma_H_rlu * randn(rng),
            Float64(Ktrue) + resolution.sigma_K_rlu * randn(rng),
            Float64(Ltrue) + resolution.sigma_L_rlu * randn(rng),
        )
    elseif resolution.momentum_mode == :cartesian_Ainv
        qtrue = rlu_to_cart(Htrue, Ktrue, Ltrue, lattice)
        qmeas = qtrue .+ resolution.sigma_Q_Ainv .* randn(rng, 3)
        Hm, Km, Lm = cart_to_rlu(qmeas, lattice)
        return (Hm, Km, Lm)
    else
        error("Unknown momentum_mode = $(resolution.momentum_mode).")
    end
end

function expand_range_by_sigma(range::Tuple{Float64, Float64}, sigma::Real, nsigma::Real)
    margin = max(0.0, Float64(nsigma) * Float64(sigma))
    return (range[1] - margin, range[2] + margin)
end

function expanded_hkl_ranges(cut::CutSpec1D, resolution::ResolutionParams, lattice::LatticeParams)
    sigmaH, sigmaK, sigmaL = momentum_sigmas_hkl(resolution, lattice)
    nsig = resolution.momentum_enabled ? resolution.sample_margin_nsigma : 0.0
    return (
        H_range = expand_range_by_sigma(cut.H_range, sigmaH, nsig),
        K_range = expand_range_by_sigma(cut.K_range, sigmaK, nsig),
        L_range = expand_range_by_sigma(cut.L_range, sigmaL, nsig),
    )
end

function uvw_sigmas_from_hkl(resolution::ResolutionParams, spec::PathCutSpec2D, lattice::LatticeParams)
    sigmaH, sigmaK, sigmaL = momentum_sigmas_hkl(resolution, lattice)
    A = projection_matrix(spec) \ Matrix{Float64}(I, 3, 3)
    Σhkl = Diagonal([sigmaH^2, sigmaK^2, sigmaL^2])
    Σuvw = A * Σhkl * transpose(A)
    return (sqrt(max(Σuvw[1,1], 0.0)), sqrt(max(Σuvw[2,2], 0.0)), sqrt(max(Σuvw[3,3], 0.0)))
end

function expanded_uvw_ranges(spec::PathCutSpec2D, resolution::ResolutionParams, lattice::LatticeParams)
    sigma_u, sigma_v, sigma_w = uvw_sigmas_from_hkl(resolution, spec, lattice)
    nsig = resolution.momentum_enabled ? resolution.sample_margin_nsigma : 0.0
    return (
        u_range = expand_range_by_sigma(spec.u_range, sigma_u, nsig),
        v_range = expand_range_by_sigma(spec.v_range, sigma_v, nsig),
        w_range = expand_range_by_sigma(spec.w_range, sigma_w, nsig),
    )
end

# -----------------------------
# Disorder sampling
# -----------------------------

function sample_disordered_parameters(model::ModelParams, disorder::DisorderParams, rng::AbstractRNG)
    if !disorder.enabled
        return model.gzz, model.J1_meV, model.J2_meV
    end

    sigmaJ1_abs = if disorder.J_units == :fractional
        disorder.sigma_J1 * abs(model.J1_meV)
    elseif disorder.J_units == :absolute_meV
        disorder.sigma_J1
    else
        error("disorder.J_units must be :fractional or :absolute_meV")
    end

    sigmaJ2_abs = if disorder.J_units == :fractional
        disorder.sigma_J2 * abs(model.J2_meV)
    elseif disorder.J_units == :absolute_meV
        disorder.sigma_J2
    else
        error("disorder.J_units must be :fractional or :absolute_meV")
    end

    zg = randn(rng)
    z1 = randn(rng)
    z2 = disorder.correlate_J1_J2 ? z1 : randn(rng)

    gzz = model.gzz + disorder.sigma_gzz * zg
    J1 = model.J1_meV + sigmaJ1_abs * z1
    J2 = model.J2_meV + sigmaJ2_abs * z2

    return gzz, J1, J2
end

# -----------------------------
# 1D simulation
# -----------------------------

function suggested_nsamples_1d(cut::CutSpec1D, lattice::LatticeParams, resolution::ResolutionParams;
                               oversample::Real = 2.0,
                               min_samples::Int = 50_000,
                               max_samples::Int = 300_000)
    sigmaH, sigmaK, sigmaL = momentum_sigmas_hkl(resolution, lattice)
    ranges = expanded_hkl_ranges(cut, resolution, lattice)

    dH = ranges.H_range[2] - ranges.H_range[1]
    dK = ranges.K_range[2] - ranges.K_range[1]
    dL = ranges.L_range[2] - ranges.L_range[1]

    stepH = sigmaH > 0 ? sigmaH / oversample : max(dH / 25.0, 1e-6)
    stepK = sigmaK > 0 ? sigmaK / oversample : max(dK / 25.0, 1e-6)
    stepL = sigmaL > 0 ? sigmaL / oversample : max(dL / 25.0, 1e-6)

    nH = max(3, ceil(Int, dH / stepH))
    nK = max(3, ceil(Int, dK / stepK))
    nL = max(3, ceil(Int, dL / stepL))

    return clamp(nH * nK * nL, min_samples, max_samples)
end

function simulate_cut_1d(cut::CutSpec1D, model::ModelParams;
                         lattice::LatticeParams = LatticeParams(),
                         disorder::DisorderParams = DisorderParams(),
                         resolution::ResolutionParams = ResolutionParams(),
                         rng::AbstractRNG = MersenneTwister(1234),
                         n_samples::Union{Nothing, Int} = nothing,
                         use_form_factor::Bool = true,
                         include_j2_formfactor::Bool = true,
                         polarization::Symbol = :transverse_c,
                         include_kfki::Bool = true,
                         intensity_scale::Real = 1.0,
                         normalize_by_samples::Bool = true)

    edges = collect(cut.E_range[1]:cut.dE_meV:cut.E_range[2])
    centers = bin_centers(edges)
    hist = zeros(Float64, length(centers))

    nsamp = isnothing(n_samples) ? suggested_nsamples_1d(cut, lattice, resolution) : n_samples
    true_ranges = expanded_hkl_ranges(cut, resolution, lattice)

    naccepted = 0
    for _ in 1:nsamp
        Htrue = rand_uniform(rng, true_ranges.H_range)
        Ktrue = rand_uniform(rng, true_ranges.K_range)
        Ltrue = rand_uniform(rng, true_ranges.L_range)

        Hm, Km, Lm = smear_hkl_to_measured(Htrue, Ktrue, Ltrue, resolution, lattice, rng)
        if !(in_range(Hm, cut.H_range) && in_range(Km, cut.K_range) && in_range(Lm, cut.L_range))
            continue
        end
        naccepted += 1

        gzz, J1, J2 = sample_disordered_parameters(model, disorder, rng)
        Etrue = dispersion_meV(Htrue, Ktrue, model; gzz = gzz, J1_meV = J1, J2_meV = J2)
        if !isfinite(Etrue) || Etrue < 0.0 || Etrue >= resolution.Ei_meV
            continue
        end

        qtrue = rlu_to_cart(Htrue, Ktrue, Ltrue, lattice)
        Qmag = norm(qtrue)
        ff = use_form_factor ? yb3_form_factor(Qmag; include_j2 = include_j2_formfactor) : 1.0
        pol = polarization_factor(qtrue; mode = polarization)
        kin = include_kfki ? kf_over_ki(Etrue, resolution.Ei_meV) : 1.0
        weight = Float64(intensity_scale) * ff^2 * pol * kin

        if resolution.energy_enabled
            sigmaE = energy_resolution_sigma_meV(Etrue, resolution.Ei_meV)
            deposit_gaussian_1d!(hist, edges, Etrue, sigmaE, weight)
        else
            deposit_delta_1d!(hist, edges, Etrue, weight)
        end
    end

    if normalize_by_samples && nsamp > 0
        hist ./= nsamp
    end

    return (
        name = cut.name,
        cut = cut,
        model = model,
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        n_samples = nsamp,
        naccepted = naccepted,
        acceptance_fraction = nsamp > 0 ? naccepted / nsamp : NaN,
        true_q_ranges = true_ranges,
        E_edges_meV = edges,
        E_centers_meV = centers,
        intensity = hist,
    )
end

function simulate_all_cuts_1d(model::ModelParams;
                               cuts::Vector{CutSpec1D} = default_cuts_1d(),
                               kwargs...)
    out = Dict{String, Any}()
    for cut in cuts
        out[cut.name] = simulate_cut_1d(cut, model; kwargs...)
    end
    return out
end

function sample_true_from_measured_hkl_1d(Hm::Real, Km::Real, Lm::Real,
                                          resolution::ResolutionParams,
                                          lattice::LatticeParams,
                                          rng::AbstractRNG)
    # Same conditional-resolution logic as the stratified 2D map.
    #
    # We sample measured Q directly inside the experimental cut volume, then
    # sample true Q from the Gaussian momentum/mosaic resolution around that
    # measured point. This includes bleed-in from outside the nominal cut without
    # the low-acceptance true-Q accept/reject step.
    if !resolution.momentum_enabled
        return (Float64(Hm), Float64(Km), Float64(Lm))
    end

    if resolution.momentum_mode == :rlu_diagonal
        return (
            Float64(Hm) + resolution.sigma_H_rlu * randn(rng),
            Float64(Km) + resolution.sigma_K_rlu * randn(rng),
            Float64(Lm) + resolution.sigma_L_rlu * randn(rng),
        )
    elseif resolution.momentum_mode == :cartesian_Ainv
        qmeas = rlu_to_cart(Hm, Km, Lm, lattice)
        qtrue = qmeas .+ resolution.sigma_Q_Ainv .* randn(rng, 3)
        Htrue, Ktrue, Ltrue = cart_to_rlu(qtrue, lattice)
        return (Htrue, Ktrue, Ltrue)
    else
        error("Unknown momentum_mode = $(resolution.momentum_mode).")
    end
end

function simulate_cut_1d_stratified(cut::CutSpec1D, model::ModelParams;
                                    lattice::LatticeParams = LatticeParams(),
                                    disorder::DisorderParams = DisorderParams(),
                                    resolution::ResolutionParams = ResolutionParams(),
                                    rng::AbstractRNG = MersenneTwister(1234),
                                    n_samples::Int = 200_000,
                                    use_form_factor::Bool = true,
                                    include_j2_formfactor::Bool = true,
                                    polarization::Symbol = :transverse_c,
                                    include_kfki::Bool = true,
                                    intensity_scale::Real = 1.0,
                                    normalize_by_samples::Bool = true)

    edges = collect(cut.E_range[1]:cut.dE_meV:cut.E_range[2])
    centers = bin_centers(edges)
    hist = zeros(Float64, length(centers))

    naccepted = 0

    for _ in 1:n_samples
        # Sample measured Q directly inside the displayed/integrated cut.
        Hm = rand_uniform(rng, cut.H_range)
        Km = rand_uniform(rng, cut.K_range)
        Lm = rand_uniform(rng, cut.L_range)

        # Draw true Q from the momentum/mosaic resolution kernel around the
        # measured Q. This is the efficient conditional form of the convolution.
        Htrue, Ktrue, Ltrue = sample_true_from_measured_hkl_1d(Hm, Km, Lm, resolution, lattice, rng)

        gzz, J1, J2 = sample_disordered_parameters(model, disorder, rng)
        Etrue = dispersion_meV(Htrue, Ktrue, model; gzz = gzz, J1_meV = J1, J2_meV = J2)

        if !isfinite(Etrue) || Etrue < 0.0 || Etrue >= resolution.Ei_meV
            continue
        end

        naccepted += 1

        qtrue = rlu_to_cart(Htrue, Ktrue, Ltrue, lattice)
        Qmag = norm(qtrue)
        ff = use_form_factor ? yb3_form_factor(Qmag; include_j2 = include_j2_formfactor) : 1.0
        pol = polarization_factor(qtrue; mode = polarization)
        kin = include_kfki ? kf_over_ki(Etrue, resolution.Ei_meV) : 1.0
        weight = Float64(intensity_scale) * ff^2 * pol * kin

        if resolution.energy_enabled
            sigmaE = energy_resolution_sigma_meV(Etrue, resolution.Ei_meV)
            deposit_gaussian_1d!(hist, edges, Etrue, sigmaE, weight)
        else
            deposit_delta_1d!(hist, edges, Etrue, weight)
        end
    end

    if normalize_by_samples && n_samples > 0
        hist ./= n_samples
    end

    return (
        name = cut.name,
        cut = cut,
        model = model,
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        n_samples = n_samples,
        naccepted = naccepted,
        acceptance_fraction = n_samples > 0 ? naccepted / n_samples : NaN,
        E_edges_meV = edges,
        E_centers_meV = centers,
        intensity = hist,
        sampling_mode = :stratified_measured_Q,
    )
end

function simulate_all_cuts_1d_stratified(model::ModelParams;
                                         cuts::Vector{CutSpec1D} = default_cuts_1d(),
                                         n_samples_per_cut::Int = 200_000,
                                         kwargs...)
    out = Dict{String, Any}()
    for cut in cuts
        out[cut.name] = simulate_cut_1d_stratified(
            cut, model;
            n_samples = n_samples_per_cut,
            kwargs...
        )
    end
    return out
end

# -----------------------------
# 2D simulation: E vs momentum along a leg
# -----------------------------

function suggested_nsamples_2d(spec::PathCutSpec2D, lattice::LatticeParams, resolution::ResolutionParams;
                               oversample::Real = 2.0,
                               min_samples::Int = 150_000,
                               max_samples::Int = 600_000)
    ranges = expanded_uvw_ranges(spec, resolution, lattice)
    sigma_u, sigma_v, sigma_w = uvw_sigmas_from_hkl(resolution, spec, lattice)

    du = ranges.u_range[2] - ranges.u_range[1]
    dv = ranges.v_range[2] - ranges.v_range[1]
    dw = ranges.w_range[2] - ranges.w_range[1]

    stepu = sigma_u > 0 ? sigma_u / oversample : max(spec.du / 2.0, 1e-6)
    stepv = sigma_v > 0 ? sigma_v / oversample : max(dv / 25.0, 1e-6)
    stepw = sigma_w > 0 ? sigma_w / oversample : max(dw / 25.0, 1e-6)

    nu = max(5, ceil(Int, du / stepu))
    nv = max(3, ceil(Int, dv / stepv))
    nw = max(3, ceil(Int, dw / stepw))

    return clamp(nu * nv * nw, min_samples, max_samples)
end

function simulate_path_map_2d(spec::PathCutSpec2D, model::ModelParams;
                              lattice::LatticeParams = LatticeParams(),
                              disorder::DisorderParams = DisorderParams(),
                              resolution::ResolutionParams = ResolutionParams(),
                              rng::AbstractRNG = MersenneTwister(1234),
                              n_samples::Union{Nothing, Int} = nothing,
                              use_form_factor::Bool = true,
                              include_j2_formfactor::Bool = true,
                              polarization::Symbol = :transverse_c,
                              include_kfki::Bool = true,
                              intensity_scale::Real = 1.0,
                              normalize_by_samples::Bool = true)

    uedges = collect(spec.u_range[1]:spec.du:spec.u_range[2])
    eedges = collect(spec.E_range[1]:spec.dE_meV:spec.E_range[2])
    ucenters = bin_centers(uedges)
    ecenters = bin_centers(eedges)
    hist = zeros(Float64, length(ucenters), length(ecenters))

    nsamp = isnothing(n_samples) ? suggested_nsamples_2d(spec, lattice, resolution) : n_samples
    true_ranges = expanded_uvw_ranges(spec, resolution, lattice)

    naccepted = 0
    for _ in 1:nsamp
        utrue = rand_uniform(rng, true_ranges.u_range)
        vtrue = rand_uniform(rng, true_ranges.v_range)
        wtrue = rand_uniform(rng, true_ranges.w_range)

        Htrue, Ktrue, Ltrue = uvw_to_hkl(utrue, vtrue, wtrue, spec)
        Hm, Km, Lm = smear_hkl_to_measured(Htrue, Ktrue, Ltrue, resolution, lattice, rng)
        umeas, vmeas, wmeas = hkl_to_uvw(Hm, Km, Lm, spec)

        if !(in_range(umeas, spec.u_range) && in_range(vmeas, spec.v_range) && in_range(wmeas, spec.w_range))
            continue
        end
        naccepted += 1

        gzz, J1, J2 = sample_disordered_parameters(model, disorder, rng)
        Etrue = dispersion_meV(Htrue, Ktrue, model; gzz = gzz, J1_meV = J1, J2_meV = J2)
        if !isfinite(Etrue) || Etrue < 0.0 || Etrue >= resolution.Ei_meV
            continue
        end

        qtrue = rlu_to_cart(Htrue, Ktrue, Ltrue, lattice)
        Qmag = norm(qtrue)
        ff = use_form_factor ? yb3_form_factor(Qmag; include_j2 = include_j2_formfactor) : 1.0
        pol = polarization_factor(qtrue; mode = polarization)
        kin = include_kfki ? kf_over_ki(Etrue, resolution.Ei_meV) : 1.0
        weight = Float64(intensity_scale) * ff^2 * pol * kin

        sigmaE = resolution.energy_enabled ? energy_resolution_sigma_meV(Etrue, resolution.Ei_meV) : 0.0
        deposit_2d_energyline!(hist, uedges, eedges, umeas, Etrue, sigmaE, weight)
    end

    if normalize_by_samples && nsamp > 0
        hist ./= nsamp
    end

    return (
        name = spec.name,
        spec = spec,
        model = model,
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        n_samples = nsamp,
        naccepted = naccepted,
        acceptance_fraction = nsamp > 0 ? naccepted / nsamp : NaN,
        true_uvw_ranges = true_ranges,
        u_edges = uedges,
        u_centers = ucenters,
        E_edges_meV = eedges,
        E_centers_meV = ecenters,
        intensity = hist,
    )
end


function sample_true_from_measured_hkl(Hm::Real, Km::Real, Lm::Real,
                                       resolution::ResolutionParams,
                                       lattice::LatticeParams,
                                       rng::AbstractRNG)
    # Conditional-resolution sampler.
    #
    # Instead of sampling true Q from a large box and rejecting most events, this
    # samples a measured Q inside the displayed histogram bin and then samples a
    # true Q from the Gaussian resolution function around that measured Q.
    #
    # For a symmetric Gaussian resolution and a locally uniform prior in Q, this
    # is equivalent to the forward true->measured convolution but is vastly more
    # efficient and gives uniform statistics in every displayed u bin.
    if !resolution.momentum_enabled
        return (Float64(Hm), Float64(Km), Float64(Lm))
    end

    if resolution.momentum_mode == :rlu_diagonal
        return (
            Float64(Hm) + resolution.sigma_H_rlu * randn(rng),
            Float64(Km) + resolution.sigma_K_rlu * randn(rng),
            Float64(Lm) + resolution.sigma_L_rlu * randn(rng),
        )
    elseif resolution.momentum_mode == :cartesian_Ainv
        qmeas = rlu_to_cart(Hm, Km, Lm, lattice)
        qtrue = qmeas .+ resolution.sigma_Q_Ainv .* randn(rng, 3)
        Htrue, Ktrue, Ltrue = cart_to_rlu(qtrue, lattice)
        return (Htrue, Ktrue, Ltrue)
    else
        error("Unknown momentum_mode = $(resolution.momentum_mode).")
    end
end

function simulate_path_map_2d_stratified(spec::PathCutSpec2D, model::ModelParams;
                                         lattice::LatticeParams = LatticeParams(),
                                         disorder::DisorderParams = DisorderParams(),
                                         resolution::ResolutionParams = ResolutionParams(),
                                         rng::AbstractRNG = MersenneTwister(1234),
                                         samples_per_u_bin::Int = 20_000,
                                         use_form_factor::Bool = true,
                                         include_j2_formfactor::Bool = true,
                                         polarization::Symbol = :transverse_c,
                                         include_kfki::Bool = true,
                                         intensity_scale::Real = 1.0,
                                         normalize_by_column_samples::Bool = true)

    uedges = collect(spec.u_range[1]:spec.du:spec.u_range[2])
    eedges = collect(spec.E_range[1]:spec.dE_meV:spec.E_range[2])
    ucenters = bin_centers(uedges)
    ecenters = bin_centers(eedges)
    hist = zeros(Float64, length(ucenters), length(ecenters))
    accepted_by_u = zeros(Int, length(ucenters))

    # No rejection is needed for u/v/w acceptance: we sample measured Q directly
    # inside each displayed bin and integration range. Resolution bleed-in is
    # included because true Q is sampled around measured Q.
    for iu in eachindex(ucenters)
        u_bin = (uedges[iu], uedges[iu + 1])

        for _ in 1:samples_per_u_bin
            umeas = rand_uniform(rng, u_bin)
            vmeas = rand_uniform(rng, spec.v_range)
            wmeas = rand_uniform(rng, spec.w_range)

            Hm, Km, Lm = uvw_to_hkl(umeas, vmeas, wmeas, spec)
            Htrue, Ktrue, Ltrue = sample_true_from_measured_hkl(Hm, Km, Lm, resolution, lattice, rng)

            gzz, J1, J2 = sample_disordered_parameters(model, disorder, rng)
            Etrue = dispersion_meV(Htrue, Ktrue, model; gzz = gzz, J1_meV = J1, J2_meV = J2)

            if !isfinite(Etrue) || Etrue < 0.0 || Etrue >= resolution.Ei_meV
                continue
            end

            accepted_by_u[iu] += 1

            qtrue = rlu_to_cart(Htrue, Ktrue, Ltrue, lattice)
            Qmag = norm(qtrue)
            ff = use_form_factor ? yb3_form_factor(Qmag; include_j2 = include_j2_formfactor) : 1.0
            pol = polarization_factor(qtrue; mode = polarization)
            kin = include_kfki ? kf_over_ki(Etrue, resolution.Ei_meV) : 1.0
            weight = Float64(intensity_scale) * ff^2 * pol * kin

            sigmaE = resolution.energy_enabled ? energy_resolution_sigma_meV(Etrue, resolution.Ei_meV) : 0.0

            # Deposit into this column only. We already conditioned on the
            # displayed measured-u bin, so no horizontal Poisson noise is added.
            if sigmaE <= 0.0 || !isfinite(sigmaE)
                ie = searchsortedlast(eedges, Etrue)
                if 1 <= ie <= size(hist, 2)
                    hist[iu, ie] += weight
                end
            else
                loE = Etrue - 6.0 * sigmaE
                hiE = Etrue + 6.0 * sigmaE
                j0 = max(1, searchsortedlast(eedges, loE))
                j1 = min(size(hist, 2), searchsortedfirst(eedges, hiE))
                for j in j0:j1
                    pbin = gaussian_bin_probability(eedges[j], eedges[j + 1], Etrue, sigmaE)
                    hist[iu, j] += weight * pbin
                end
            end
        end

        if normalize_by_column_samples && samples_per_u_bin > 0
            hist[iu, :] ./= samples_per_u_bin
        end
    end

    n_samples = samples_per_u_bin * length(ucenters)
    naccepted = sum(accepted_by_u)

    return (
        name = spec.name,
        spec = spec,
        model = model,
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        n_samples = n_samples,
        naccepted = naccepted,
        acceptance_fraction = n_samples > 0 ? naccepted / n_samples : NaN,
        samples_per_u_bin = samples_per_u_bin,
        accepted_by_u = accepted_by_u,
        u_edges = uedges,
        u_centers = ucenters,
        E_edges_meV = eedges,
        E_centers_meV = ecenters,
        intensity = hist,
    )
end

# -----------------------------
# CSV export
# -----------------------------

function write_cut_csv(filename::AbstractString, result)
    open(filename, "w") do io
        println(io, "E_meV,intensity")
        for (E, I) in zip(result.E_centers_meV, result.intensity)
            println(io, @sprintf("%.8g,%.12g", E, I))
        end
    end
    return filename
end

function write_all_csv_1d(prefix::AbstractString, results::Dict{String, Any})
    files = String[]
    for (name, result) in results
        safe_name = replace(name, r"[^A-Za-z0-9_\-]" => "_")
        push!(files, write_cut_csv("$(prefix)_$(safe_name).csv", result))
    end
    return files
end

function write_map_csv_2d(filename::AbstractString, result)
    open(filename, "w") do io
        println(io, join(vcat("u_over_E", string.(result.E_centers_meV)), ","))
        for i in eachindex(result.u_centers)
            row = [@sprintf("%.8g", result.u_centers[i])]
            append!(row, [@sprintf("%.12g", result.intensity[i, j]) for j in eachindex(result.E_centers_meV)])
            println(io, join(row, ","))
        end
    end
    return filename
end

# -----------------------------
# Print helpers
# -----------------------------

function print_model(model::ModelParams)
    println("ModelParams")
    println("-----------")
    println(@sprintf("B     = %.6g T", model.B_T))
    println(@sprintf("gzz   = %.6g", model.gzz))
    println(@sprintf("J1    = %.6g meV", model.J1_meV))
    println(@sprintf("J2    = %.6g meV", model.J2_meV))
    println(@sprintf("J2/J1 = %.6g", model.J2_meV / model.J1_meV))
    println(@sprintf("S     = %.6g", model.S))
end

function high_symmetry_energies(model::ModelParams)
    EGamma = dispersion_meV(0.0, 0.0, model)
    EK = dispersion_meV(1.0 / 3.0, 1.0 / 3.0, model)
    EM = dispersion_meV(0.5, 0.0, model)
    return (EGamma_meV = EGamma, EK_meV = EK, EM_meV = EM)
end

function print_momentum_resolution(resolution::ResolutionParams, lattice::LatticeParams = LatticeParams())
    println("Momentum resolution / mosaic")
    println("----------------------------")
    println("enabled = $(resolution.momentum_enabled)")
    println("mode    = $(resolution.momentum_mode)")
    println()

    if resolution.momentum_mode == :rlu_diagonal
        println(@sprintf("sigma_H = %.6g rlu, FWHM_H = %.6g rlu", resolution.sigma_H_rlu, resolution.sigma_H_rlu * SIGMA_TO_FWHM))
        println(@sprintf("sigma_K = %.6g rlu, FWHM_K = %.6g rlu", resolution.sigma_K_rlu, resolution.sigma_K_rlu * SIGMA_TO_FWHM))
        println(@sprintf("sigma_L = %.6g rlu, FWHM_L = %.6g rlu", resolution.sigma_L_rlu, resolution.sigma_L_rlu * SIGMA_TO_FWHM))
    else
        Bmat = rlu_basis_matrix(lattice)
        sigmaH, sigmaK, sigmaL = momentum_sigmas_hkl(resolution, lattice)
        println(@sprintf("sigma_Q = %.6g Å^-1", resolution.sigma_Q_Ainv))
        println(@sprintf("|a*|    = %.6g Å^-1", norm(Bmat[:, 1])))
        println(@sprintf("|b*|    = %.6g Å^-1", norm(Bmat[:, 2])))
        println(@sprintf("|c*|    = %.6g Å^-1", norm(Bmat[:, 3])))
        println(@sprintf("sigma_H = %.6g rlu, FWHM_H = %.6g rlu", sigmaH, sigmaH * SIGMA_TO_FWHM))
        println(@sprintf("sigma_K = %.6g rlu, FWHM_K = %.6g rlu", sigmaK, sigmaK * SIGMA_TO_FWHM))
        println(@sprintf("sigma_L = %.6g rlu, FWHM_L = %.6g rlu", sigmaL, sigmaL * SIGMA_TO_FWHM))
    end

    println()
    println(@sprintf("sample margin = %.3g sigma outside each experimental bound", resolution.sample_margin_nsigma))
end

function ordered_result_names_1d(results::Dict{String, Any})
    preferred = [
        "Gamma_cut_center_0_1_0",
        "M_label_cut_center_1over3_1over3_0",
        "K_label_cut_center_1over2_0_0",
    ]
    names = String[]
    for name in preferred
        if haskey(results, name)
            push!(names, name)
        end
    end
    for name in sort!(collect(keys(results)))
        if !(name in names)
            push!(names, name)
        end
    end
    return names
end

function print_cut_sampling_summary(results::Dict{String, Any})
    println("Sampling summary (1D cuts)")
    println("--------------------------")
    for name in ordered_result_names_1d(results)
        r = results[name]
        println(name)
        println(@sprintf("  n_samples           = %d", r.n_samples))
        println(@sprintf("  naccepted           = %d", r.naccepted))
        println(@sprintf("  acceptance_fraction = %.6g", r.acceptance_fraction))
        println(@sprintf("  true H range        = [%.6g, %.6g]", r.true_q_ranges.H_range...))
        println(@sprintf("  true K range        = [%.6g, %.6g]", r.true_q_ranges.K_range...))
        println(@sprintf("  true L range        = [%.6g, %.6g]", r.true_q_ranges.L_range...))
    end
end

# -----------------------------
# GLMakie plotting: 1D
# -----------------------------

function cut_pretty_title_1d(name::AbstractString)
    if name == "Gamma_cut_center_0_1_0"
        return "Gamma cut  (center ~ (0,1,0))"
    elseif name == "M_label_cut_center_1over3_1over3_0"
        return "Experimental 'M' cut  (analytic K ~ (1/3,1/3,0))"
    elseif name == "K_label_cut_center_1over2_0_0"
        return "Experimental 'K' cut  (analytic M ~ (1/2,0,0))"
    else
        return String(name)
    end
end

function plot_results_1d(results::Dict{String, Any}; figure_title::AbstractString = "YZGO analytic cut simulation")
    names = ordered_result_names_1d(results)
    n = length(names)
    fig = Figure(size = (950, 270 * max(n, 1)))
    Label(fig[0, 1], figure_title; fontsize = 22)

    axes = Axis[]
    for (i, name) in enumerate(names)
        result = results[name]
        ax = Axis(fig[i, 1], xlabel = "Energy transfer (meV)", ylabel = "Intensity (arb. units)", title = cut_pretty_title_1d(name))
        lines!(ax, result.E_centers_meV, result.intensity, linewidth = 2.0)
        xlims!(ax, result.cut.E_range...)
        push!(axes, ax)
    end
    linkxaxes!(axes...)
    display(fig)
    return fig
end

function compare_results_plot_1d(result_sets::Vector{Pair{String, Dict{String, Any}}};
                                 figure_title::AbstractString = "YZGO analytic cut simulation")
    isempty(result_sets) && error("result_sets cannot be empty")
    names = ordered_result_names_1d(result_sets[1].second)
    n = length(names)

    fig = Figure(size = (980, 280 * max(n, 1)))
    Label(fig[0, 1], figure_title; fontsize = 22)
    axes = Axis[]
    for (i, cutname) in enumerate(names)
        ax = Axis(fig[i, 1], xlabel = "Energy transfer (meV)", ylabel = "Intensity (arb. units)", title = cut_pretty_title_1d(cutname))
        push!(axes, ax)
        for (setlabel, results) in result_sets
            if haskey(results, cutname)
                result = results[cutname]
                lines!(ax, result.E_centers_meV, result.intensity, linewidth = 2.0, label = setlabel)
            end
        end
        axislegend(ax; position = :rt)
    end
    linkxaxes!(axes...)
    display(fig)
    return fig
end

# -----------------------------
# GLMakie plotting: 2D path maps
# -----------------------------

function plot_path_map_comparison_2d(clean_result, disorder_result;
                                     figure_title::AbstractString = "YZGO analytic E vs momentum map")
    fig = Figure(size = (1350, 540))
    Label(fig[0, 1:2], figure_title; fontsize = 22)

    xtickvals = [-1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0]
    xticklabels = ["K₁", "M₁", "K", "Γ₁"]
    colmax = maximum(vcat(vec(clean_result.intensity), vec(disorder_result.intensity)))
    crange = (0.0, colmax > 0 ? colmax : 1.0)

    map_results = [(clean_result, "Without disorder"), (disorder_result, "With disorder")]
    for (icol, (res, title)) in enumerate(map_results)
        ax = Axis(fig[1, icol],
            xlabel = "u along K₁ → M₁ → K → Γ₁",
            ylabel = "Energy transfer (meV)",
            title = title,
            xticks = (xtickvals, xticklabels),
        )

        hm = heatmap!(ax, res.u_centers, res.E_centers_meV, res.intensity; colorrange = crange)
        vlines!(ax, xtickvals, color = (:white, 0.8), linestyle = :dash, linewidth = 1.0)
        xlims!(ax, res.spec.u_range...)
        ylims!(ax, res.spec.E_range...)
        Colorbar(fig[1, icol + 2], hm, label = "Intensity (arb. units)")
    end

    display(fig)
    return fig
end

# -----------------------------
# Demo / convenience runners
# -----------------------------

function demo_defaults()
    model = fit_model_from_GKM(
        1.825,  # E_Gamma at 9 T
        1.025,  # E at analytic K = (1/3,1/3)
        1.075,  # E at analytic M = (1/2,0)
        9.0
    )

    lattice = LatticeParams(a_A = 3.376, c_A = 21.96)

    disorder = DisorderParams(
        enabled = true,
        sigma_gzz = 0.25, #initial parameter for gzz dispersive
        sigma_J1 = 1.0 / 3.0,
        sigma_J2 = 1.0 / 3.0,
        J_units = :fractional,
        correlate_J1_J2 = false,
    )

    disorder_off = DisorderParams(enabled = false)

    resolution = ResolutionParams(
        Ei_meV = 4.65,
        energy_enabled = true,
        momentum_enabled = true,
        momentum_mode = :rlu_diagonal,
        sigma_H_rlu = 0.037,
        sigma_K_rlu = 0.037,
        sigma_L_rlu = 0.252,
        sample_margin_nsigma = 3.0,
    )

    return (; model, lattice, disorder, disorder_off, resolution)
end

function run_demo_1d()
    pars = demo_defaults()
    model, lattice, disorder, resolution = pars.model, pars.lattice, pars.disorder, pars.resolution

    print_model(model)
    println()
    println("High-symmetry energies from model:")
    println(high_symmetry_energies(model))
    println()
    print_momentum_resolution(resolution, lattice)
    println()

    results = simulate_all_cuts_1d(
        model;
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        rng = MersenneTwister(2026),
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )

    print_cut_sampling_summary(results)
    println()

    files = write_all_csv_1d("YZGO_analytic_1d", results)
    println("Wrote CSV files:")
    foreach(println, files)

    fig = plot_results_1d(results; figure_title = "YZGO analytic 1D cuts")
    return (; results = results, figure = fig)
end

function run_demo_2d()
    pars = demo_defaults()
    model, lattice = pars.model, pars.lattice
    disorder, disorder_off, resolution = pars.disorder, pars.disorder_off, pars.resolution
    spec = default_leg_cut_2d()

    print_model(model)
    println()
    println("High-symmetry energies from model:")
    println(high_symmetry_energies(model))
    println()
    print_momentum_resolution(resolution, lattice)
    println()
    println("2D path spec:")
    println(spec)
    println()

    clean_map = simulate_path_map_2d(
        spec, model;
        lattice = lattice,
        disorder = disorder_off,
        resolution = resolution,
        rng = MersenneTwister(2026),
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )

    disorder_map = simulate_path_map_2d(
        spec, model;
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        rng = MersenneTwister(2027),
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )

    println(@sprintf("2D clean:    n_samples = %d, naccepted = %d, acceptance_fraction = %.6g", clean_map.n_samples, clean_map.naccepted, clean_map.acceptance_fraction))
    println(@sprintf("2D disorder: n_samples = %d, naccepted = %d, acceptance_fraction = %.6g", disorder_map.n_samples, disorder_map.naccepted, disorder_map.acceptance_fraction))

    file_clean = write_map_csv_2d("YZGO_analytic_2d_clean.csv", clean_map)
    file_dis = write_map_csv_2d("YZGO_analytic_2d_disorder.csv", disorder_map)
    println("Wrote CSV files:")
    println(file_clean)
    println(file_dis)

    fig = plot_path_map_comparison_2d(clean_map, disorder_map; figure_title = "YZGO analytic 2D map: K₁ → M₁ → K → Γ₁")
    return (; clean = clean_map, disorder = disorder_map, figure = fig)
end

function run_demo_1d_stratified(; n_samples_per_cut::Int = 500_000)
    pars = demo_defaults()
    model, lattice, disorder, resolution = pars.model, pars.lattice, pars.disorder, pars.resolution

    print_model(model)
    println()
    println("High-symmetry energies from model:")
    println(high_symmetry_energies(model))
    println()
    print_momentum_resolution(resolution, lattice)
    println()
    println(@sprintf("Stratified 1D sampling with %d samples per cut", n_samples_per_cut))
    println()

    results = simulate_all_cuts_1d_stratified(
        model;
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        rng = MersenneTwister(2026),
        n_samples_per_cut = n_samples_per_cut,
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )

    println("Sampling summary (1D stratified)")
    println("--------------------------------")
    for name in ordered_result_names_1d(results)
        r = results[name]
        println(name)
        println(@sprintf("  n_samples           = %d", r.n_samples))
        println(@sprintf("  naccepted           = %d", r.naccepted))
        println(@sprintf("  acceptance_fraction = %.6g", r.acceptance_fraction))
    end
    println()

    files = write_all_csv_1d("YZGO_analytic_1d_stratified", results)
    println("Wrote CSV files:")
    foreach(println, files)

    fig = plot_results_1d(results; figure_title = "YZGO analytic 1D cuts, stratified")
    return (; results = results, figure = fig)
end

function run_demo_2d_stratified(; samples_per_u_bin::Int = 20_000)
    pars = demo_defaults()
    model, lattice = pars.model, pars.lattice
    disorder, disorder_off, resolution = pars.disorder, pars.disorder_off, pars.resolution
    spec = default_leg_cut_2d()

    print_model(model)
    println()
    println("High-symmetry energies from model:")
    println(high_symmetry_energies(model))
    println()
    print_momentum_resolution(resolution, lattice)
    println()
    println("2D path spec:")
    println(spec)
    println()
    println(@sprintf("Stratified 2D sampling with %d samples per u bin", samples_per_u_bin))
    println()

    clean_map = simulate_path_map_2d_stratified(
        spec, model;
        lattice = lattice,
        disorder = disorder_off,
        resolution = resolution,
        rng = MersenneTwister(2026),
        samples_per_u_bin = samples_per_u_bin,
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )

    disorder_map = simulate_path_map_2d_stratified(
        spec, model;
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        rng = MersenneTwister(2027),
        samples_per_u_bin = samples_per_u_bin,
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )

    println(@sprintf("2D clean:    n_samples = %d, naccepted = %d, acceptance_fraction = %.6g",
                     clean_map.n_samples, clean_map.naccepted, clean_map.acceptance_fraction))
    println(@sprintf("2D disorder: n_samples = %d, naccepted = %d, acceptance_fraction = %.6g",
                     disorder_map.n_samples, disorder_map.naccepted, disorder_map.acceptance_fraction))

    file_clean = write_map_csv_2d("YZGO_analytic_2d_clean_stratified.csv", clean_map)
    file_dis = write_map_csv_2d("YZGO_analytic_2d_disorder_stratified.csv", disorder_map)
    println("Wrote CSV files:")
    println(file_clean)
    println(file_dis)

    fig = plot_path_map_comparison_2d(clean_map, disorder_map;
        figure_title = "YZGO analytic 2D map, stratified: K₁ → M₁ → K → Γ₁")

    return (; clean = clean_map, disorder = disorder_map, figure = fig)
end


# =============================================================================
# Section copied from YZGO_overlay_experiment_model_1d.jl
# =============================================================================

# YZGO_overlay_experiment_model_1d.jl
#
# Small bridge script to overplot the analytical 1D model results on top of the
# experimental 1D scans loaded by yzgo_plot_1d_scans.jl.
#
# Intended workflow from Julia:
#
#   include("yzgo_plot_1d_scans.jl")
#   include("YZGO_analytic_1d2d_glmakie_stratified_1d2d.jl")
#   include("YZGO_overlay_experiment_model_1d.jl")
#
#   # Option A: make fresh stratified model cuts, then overlay them.
#   pars = demo_defaults()
#   model_results = simulate_all_cuts_1d_stratified(
#       pars.model;
#       lattice = pars.lattice,
#       disorder = pars.disorder,
#       resolution = pars.resolution,
#       rng = MersenneTwister(2026),
#       n_samples_per_cut = 500_000,
#       use_form_factor = true,
#       include_kfki = true,
#       polarization = :transverse_c,
#   )
#
#   fig, scales = plot_experiment_model_overlay_1d(
#       scans_bgsub,
#       ["model" => model_results];
#       field_T = 9.0,
#       scale_mode = :global,
#   )
#
#   # Option B: convenience wrapper that simulates the default 9 T model first.
#   out = run_overlay_demo_1d(; data = :bgsub, field_T = 9.0)
#
# Notes:
#   - Experimental data are accessed as scans_or_scans_bgsub[qtag][field_T].
#   - Model results are accessed as model_results[model_key].
#   - A multiplicative scale factor is fitted because the experimental data are
#     arbitrary units.

using GLMakie
using Printf
using Statistics

const OVERLAY_QTAG_TO_MODELKEY = Dict(
    "0_1_0"        => "Gamma_cut_center_0_1_0",
    "0p33_0p33_0"  => "M_label_cut_center_1over3_1over3_0",
    "0p5_0_0"      => "K_label_cut_center_1over2_0_0",
)

const OVERLAY_QTAG_TITLES = Dict(
    "0_1_0"        => "Q = (0, 1, 0)  /  Γ-like cut",
    "0p33_0p33_0"  => "Q = (0.33, 0.33, 0)  /  analytic K-like cut",
    "0p5_0_0"      => "Q = (0.5, 0, 0)  /  analytic M-like cut",
)

# Optional default scale-fit windows. These are intentionally broad and can be
# overridden by passing fit_windows_by_q = Dict(...). Set fit_windows_by_q=nothing
# to use the model-intensity mask only.
const DEFAULT_OVERLAY_SCALE_WINDOWS = Dict(
    "0_1_0"        => [(1.0, 2.6)],
    "0p33_0p33_0"  => [(0.4, 1.9)],
    "0p5_0_0"      => [(0.4, 1.9)],
)

function _check_overlay_prereqs()
    isdefined(Main, :scans) || @warn "`scans` is not defined. Run include(\"yzgo_plot_1d_scans.jl\") first if you want raw data."
    isdefined(Main, :scans_bgsub) || @warn "`scans_bgsub` is not defined. Run the experimental script first if you want background-subtracted data."
    return nothing
end

function _sort_xy_local(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real})
    length(xin) == length(yin) || error("x and y must have same length")
    p = sortperm(Float64.(xin))
    return Float64.(xin[p]), Float64.(yin[p])
end

function interp1_linear_local(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real}, xout::AbstractVector{<:Real})
    x, y = _sort_xy_local(xin, yin)
    n = length(x)
    n >= 2 || error("Need at least two model points for interpolation")
    any(diff(x) .<= 0) && error("Interpolation x-values must be strictly increasing")

    out = Vector{Float64}(undef, length(xout))
    for (i, xr) in enumerate(xout)
        x0 = Float64(xr)
        if x0 < x[1] || x0 > x[end]
            out[i] = NaN
            continue
        end
        j = searchsortedlast(x, x0)
        j = clamp(j, 1, n - 1)
        t = (x0 - x[j]) / (x[j + 1] - x[j])
        out[i] = (1.0 - t) * y[j] + t * y[j + 1]
    end
    return out
end

function _window_mask(E::AbstractVector{<:Real}, windows)
    mask = trues(length(E))
    if windows === nothing
        return mask
    end
    mask .= false
    for (lo, hi) in windows
        mask .|= (E .>= lo) .& (E .<= hi)
    end
    return mask
end

function _model_result_for_qtag(model_results, qtag::AbstractString;
                                qtag_to_modelkey = OVERLAY_QTAG_TO_MODELKEY)
    key = qtag_to_modelkey[String(qtag)]
    haskey(model_results, key) || error("model_results does not contain key $key for qtag $qtag")
    return model_results[key]
end

function _scan_for_qtag_field(data_scans, qtag::AbstractString, field_T::Real)
    haskey(data_scans, qtag) || error("data_scans does not contain qtag $qtag")
    byfield = data_scans[qtag]
    B = Float64(field_T)
    if haskey(byfield, B)
        return byfield[B]
    end

    # Tolerant fallback for Float64 dictionary keys.
    matches = [k for k in keys(byfield) if abs(k - B) < 1e-8]
    isempty(matches) && error("No field $(field_T) T for qtag $qtag. Available fields: $(sort(collect(keys(byfield))))")
    return byfield[first(matches)]
end

function _scale_fit_vectors(data_scans, model_results;
                            field_T::Real = 9.0,
                            qtags = OVERLAY_QTAG_TO_MODELKEY |> keys |> collect,
                            qtag_to_modelkey = OVERLAY_QTAG_TO_MODELKEY,
                            fit_windows_by_q = DEFAULT_OVERLAY_SCALE_WINDOWS,
                            min_model_fraction::Real = 0.03,
                            use_errors::Bool = true)
    y_all = Float64[]
    m_all = Float64[]
    w_all = Float64[]

    for qtag in qtags
        haskey(qtag_to_modelkey, qtag) || continue
        if !haskey(data_scans, qtag)
            @warn "Skipping qtag $qtag because it is absent from data_scans."
            continue
        end

        scan = _scan_for_qtag_field(data_scans, qtag, field_T)
        mres = _model_result_for_qtag(model_results, qtag; qtag_to_modelkey = qtag_to_modelkey)
        model_on_data = interp1_linear_local(mres.E_centers_meV, mres.intensity, scan.energy)

        finite_mask = isfinite.(scan.energy) .& isfinite.(scan.intensity) .& isfinite.(model_on_data)
        if use_errors && hasproperty(scan, :error)
            finite_mask .&= isfinite.(scan.error) .& (scan.error .> 0)
        end

        maxmodel = maximum(abs.(model_on_data[finite_mask]); init = 0.0)
        model_mask = maxmodel > 0 ? abs.(model_on_data) .>= min_model_fraction * maxmodel : trues(length(model_on_data))

        windows = fit_windows_by_q === nothing ? nothing : get(fit_windows_by_q, qtag, nothing)
        fit_mask = finite_mask .& model_mask .& _window_mask(scan.energy, windows)

        if count(fit_mask) < 3
            @warn "Very few scale-fit points for qtag $qtag; falling back to finite/model mask."
            fit_mask = finite_mask .& model_mask
        end

        if count(fit_mask) < 3
            @warn "Skipping qtag $qtag in scale fit; insufficient usable points."
            continue
        end

        append!(y_all, Float64.(scan.intensity[fit_mask]))
        append!(m_all, Float64.(model_on_data[fit_mask]))

        if use_errors && hasproperty(scan, :error)
            append!(w_all, 1.0 ./ max.(Float64.(scan.error[fit_mask]), eps(Float64)).^2)
        else
            append!(w_all, ones(count(fit_mask)))
        end
    end

    return y_all, m_all, w_all
end

function best_scale_factor(data_scans, model_results;
                           field_T::Real = 9.0,
                           qtags = collect(keys(OVERLAY_QTAG_TO_MODELKEY)),
                           qtag_to_modelkey = OVERLAY_QTAG_TO_MODELKEY,
                           fit_windows_by_q = DEFAULT_OVERLAY_SCALE_WINDOWS,
                           min_model_fraction::Real = 0.03,
                           use_errors::Bool = true,
                           positive_only::Bool = true)
    y, m, w = _scale_fit_vectors(data_scans, model_results;
        field_T = field_T,
        qtags = qtags,
        qtag_to_modelkey = qtag_to_modelkey,
        fit_windows_by_q = fit_windows_by_q,
        min_model_fraction = min_model_fraction,
        use_errors = use_errors,
    )

    length(y) >= 3 || error("Not enough points to fit a scale factor")
    denom = sum(w .* m .* m)
    denom > 0 || error("Cannot fit scale factor because model vector is zero")
    s = sum(w .* y .* m) / denom
    return positive_only ? max(s, 0.0) : s
end

function best_scale_factors_per_q(data_scans, model_results;
                                  field_T::Real = 9.0,
                                  qtags = collect(keys(OVERLAY_QTAG_TO_MODELKEY)),
                                  qtag_to_modelkey = OVERLAY_QTAG_TO_MODELKEY,
                                  fit_windows_by_q = DEFAULT_OVERLAY_SCALE_WINDOWS,
                                  min_model_fraction::Real = 0.03,
                                  use_errors::Bool = true,
                                  positive_only::Bool = true)
    out = Dict{String, Float64}()
    for qtag in qtags
        y, m, w = _scale_fit_vectors(data_scans, model_results;
            field_T = field_T,
            qtags = [qtag],
            qtag_to_modelkey = qtag_to_modelkey,
            fit_windows_by_q = fit_windows_by_q,
            min_model_fraction = min_model_fraction,
            use_errors = use_errors,
        )
        if length(y) < 3
            @warn "Not enough points to fit per-q scale for $qtag; using NaN."
            out[qtag] = NaN
            continue
        end
        denom = sum(w .* m .* m)
        s = denom > 0 ? sum(w .* y .* m) / denom : NaN
        out[qtag] = positive_only ? max(s, 0.0) : s
    end
    return out
end

function _normalize_model_sets(model_sets)
    if model_sets isa Dict
        return [String(k) => v for (k, v) in model_sets]
    elseif model_sets isa Pair
        return [String(model_sets.first) => model_sets.second]
    elseif model_sets isa AbstractVector
        return [String(p.first) => p.second for p in model_sets]
    else
        error("model_sets must be a Pair, Dict, or vector of Pairs, e.g. [\"model\" => results]")
    end
end

function _scale_lookup(scales, label::AbstractString, qtag::AbstractString)
    if scales isa Number
        return Float64(scales)
    elseif scales isa Dict
        if haskey(scales, label)
            s = scales[label]
            if s isa Number
                return Float64(s)
            elseif s isa Dict
                return Float64(s[qtag])
            else
                error("Unsupported scale entry for label $label")
            end
        elseif haskey(scales, qtag)
            return Float64(scales[qtag])
        else
            error("Scale dictionary lacks label $label or qtag $qtag")
        end
    else
        error("Unsupported scales object")
    end
end

function compute_overlay_scales(data_scans, model_sets;
                                field_T::Real = 9.0,
                                scale = :auto,
                                scale_mode::Symbol = :global,
                                qtags = collect(keys(OVERLAY_QTAG_TO_MODELKEY)),
                                qtag_to_modelkey = OVERLAY_QTAG_TO_MODELKEY,
                                fit_windows_by_q = DEFAULT_OVERLAY_SCALE_WINDOWS,
                                min_model_fraction::Real = 0.03,
                                use_errors::Bool = true)
    sets = _normalize_model_sets(model_sets)

    if scale !== :auto
        return Dict(label => scale for (label, _) in sets)
    end

    scales = Dict{String, Any}()
    for (label, results) in sets
        if scale_mode == :global
            scales[label] = best_scale_factor(data_scans, results;
                field_T = field_T,
                qtags = qtags,
                qtag_to_modelkey = qtag_to_modelkey,
                fit_windows_by_q = fit_windows_by_q,
                min_model_fraction = min_model_fraction,
                use_errors = use_errors,
            )
        elseif scale_mode == :per_q
            scales[label] = best_scale_factors_per_q(data_scans, results;
                field_T = field_T,
                qtags = qtags,
                qtag_to_modelkey = qtag_to_modelkey,
                fit_windows_by_q = fit_windows_by_q,
                min_model_fraction = min_model_fraction,
                use_errors = use_errors,
            )
        else
            error("scale_mode must be :global or :per_q")
        end
    end
    return scales
end

function print_overlay_scales(scales)
    println("Overlay scale factors")
    println("---------------------")
    for (label, s) in scales
        if s isa Number
            println(@sprintf("%s: %.8g", label, s))
        elseif s isa Dict
            println(label)
            for qtag in sort(collect(keys(s)))
                println(@sprintf("  %s: %.8g", qtag, s[qtag]))
            end
        else
            println("$label: $s")
        end
    end
end

function plot_experiment_model_overlay_1d(data_scans,
                                          model_sets;
                                          field_T::Real = 9.0,
                                          scale = :auto,
                                          scale_mode::Symbol = :global,
                                          qtags = ["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                                          qtag_to_modelkey = OVERLAY_QTAG_TO_MODELKEY,
                                          fit_windows_by_q = DEFAULT_OVERLAY_SCALE_WINDOWS,
                                          min_model_fraction::Real = 0.03,
                                          use_errors_for_scale::Bool = true,
                                          model_linewidth::Real = 3.0,
                                          data_markersize::Real = 8.0,
                                          figure_title::AbstractString = "YZGO 1D experimental data vs analytical model",
                                          save_png::Bool = false,
                                          outpath::Union{Nothing, AbstractString} = nothing)
    _check_overlay_prereqs()
    sets = _normalize_model_sets(model_sets)

    scales = compute_overlay_scales(data_scans, sets;
        field_T = field_T,
        scale = scale,
        scale_mode = scale_mode,
        qtags = qtags,
        qtag_to_modelkey = qtag_to_modelkey,
        fit_windows_by_q = fit_windows_by_q,
        min_model_fraction = min_model_fraction,
        use_errors = use_errors_for_scale,
    )

    print_overlay_scales(scales)

    fig = Figure(size = (1450, 460))
    Label(fig[0, :], figure_title * @sprintf("  (%.0f T)", field_T); fontsize = 20)

    for (icol, qtag) in enumerate(qtags)
        haskey(data_scans, qtag) || continue
        scan = _scan_for_qtag_field(data_scans, qtag, field_T)
        title = get(OVERLAY_QTAG_TITLES, qtag, qtag)

        ax = Axis(fig[1, icol];
            title = title,
            xlabel = "ΔE (meV)",
            ylabel = icol == 1 ? "Intensity (arb. units)" : "",
        )

        scatter!(ax, scan.energy, scan.intensity; label = @sprintf("data %.0f T", field_T), markersize = data_markersize)
        if hasproperty(scan, :error)
            errorbars!(ax, scan.energy, scan.intensity, scan.error; whiskerwidth = 4)
        end

        for (label, results) in sets
            mres = _model_result_for_qtag(results, qtag; qtag_to_modelkey = qtag_to_modelkey)
            s = _scale_lookup(scales, label, qtag)
            lines!(ax, mres.E_centers_meV, s .* mres.intensity; label = "model: $label", linewidth = model_linewidth)
        end

        if hasproperty(scan, :meta)
            xlims!(ax, scan.meta.de_range...)
        else
            xlims!(ax, minimum(scan.energy), maximum(scan.energy))
        end

        axislegend(ax; position = :rt, framevisible = false)
    end

    if save_png
        resolved_outpath = outpath
        if resolved_outpath === nothing
            outdir = isdefined(Main, :OUTDIR) ? Main.OUTDIR : pwd()
            resolved_outpath = joinpath(outdir, @sprintf("YZGO_1d_experiment_model_overlay_%.0fT.png", field_T))
        end
        save(resolved_outpath, fig)
        println("Saved overlay figure to: ", resolved_outpath)
    end

    display(fig)
    return fig, scales
end

function make_default_overlay_model_results(; disorder_enabled::Bool = true,
                                            n_samples_per_cut::Int = 500_000,
                                            seed::Int = 2026)
    isdefined(Main, :demo_defaults) || error("Run include(\"YZGO_analytic_1d2d_glmakie_stratified_1d2d.jl\") first.")
    pars = demo_defaults()
    disorder = disorder_enabled ? pars.disorder : pars.disorder_off
    return simulate_all_cuts_1d_stratified(
        pars.model;
        lattice = pars.lattice,
        disorder = disorder,
        resolution = pars.resolution,
        rng = MersenneTwister(seed),
        n_samples_per_cut = n_samples_per_cut,
        use_form_factor = true,
        include_kfki = true,
        polarization = :transverse_c,
    )
end

function run_overlay_demo_1d(; data::Symbol = :bgsub,
                             field_T::Real = 9.0,
                             disorder_enabled::Bool = true,
                             n_samples_per_cut::Int = 500_000,
                             scale = :auto,
                             scale_mode::Symbol = :global,
                             save_png::Bool = false)
    data_scans = if data == :bgsub
        isdefined(Main, :scans_bgsub) || error("scans_bgsub is not defined. Run include(\"yzgo_plot_1d_scans.jl\") first.")
        Main.scans_bgsub
    elseif data == :raw
        isdefined(Main, :scans) || error("scans is not defined. Run include(\"yzgo_plot_1d_scans.jl\") first.")
        Main.scans
    else
        error("data must be :bgsub or :raw")
    end

    model_results = make_default_overlay_model_results(
        disorder_enabled = disorder_enabled,
        n_samples_per_cut = n_samples_per_cut,
    )

    label = disorder_enabled ? "disorder" : "clean"
    fig, scales = plot_experiment_model_overlay_1d(
        data_scans,
        [label => model_results];
        field_T = field_T,
        scale = scale,
        scale_mode = scale_mode,
        save_png = save_png,
        figure_title = data == :bgsub ? "YZGO background-subtracted data vs analytical model" : "YZGO raw data vs analytical model",
    )

    return (; fig, scales, model_results)
end


# =============================================================================
# Monolithic 1D fitting layer: two-kernel model
# =============================================================================
#
# This section combines the data-loading/background-subtraction workflow,
# the analytical stratified measured-Q model, and a two-component fitting model.
#
# Fitted model before global scale:
#
#   model(E,Q) = I_disp(E,Q; gzz,J1,J2,sigma_gzz,sigma_J1,sigma_J2)
#              + r2 * I_flat(E,Q; gzz2,sigma_gzz2,J1=J2=0)
#
# The final plotted/fitted model is:
#
#   global_scale * model(E,Q)
#
# There is deliberately no Q=(0,1,0)-specific relative scale in this version.
#
# Default script behavior:
#   - loads the 1D scans from BASEDIR
#   - builds the same tail/background-subtracted scans as yzgo_plot_1d_scans.jl
#   - fits the 9 T high-field data only
#   - fits dispersive-kernel gzz, J1, J2, sigma_gzz, sigma_J1, sigma_J2
#   - fits non-dispersive-kernel gzz2 and sigma_gzz2
#   - fits second_kernel_relative_intensity = r2
#   - solves or fits one global scale factor
#   - masks the requested problematic-background regions
#   - plots dispersive, non-dispersive, and total model curves separately
#
# To run from Julia:
#   julia YZGO_fit_experiment_model_1d_two_kernel.jl
#
# To customize interactively:
#   include("YZGO_fit_experiment_model_1d_two_kernel.jl")
#   result = run_yzgo_fit_1d(fields_T=[9.0, 14.0], n_samples_per_cut=200_000)
#
# Notes:
#   * The objective is Monte-Carlo based, but fixed random seeds are used for
#     each evaluation so the optimizer sees a deterministic function.
#   * For faster exploratory fits, reduce n_samples_per_cut and maxiters.
#   * For final reporting, increase final_n_samples_per_cut.

using Dates

try
    @eval using Optim
catch err
    @error "This monolithic fitting script requires Optim.jl. Install with: import Pkg; Pkg.add(\"Optim\")" exception=(err, catch_backtrace())
    rethrow(err)
end

const YZGO_FIT_WINDOWS_1D = Dict(
    "0_1_0"        => [(1.1, 3.0)],
    "0p33_0p33_0" => [(0.75, 3.0)],#[(0.75, 1.6), (2.3, 3.0)],
    "0p5_0_0"     => [(0.75, 3.0)],#[(0.75, 1.6), (2.3, 3.0)],
)

const YZGO_FIT_PLOT_YLIMS = (-0.0005, 0.003)

Base.@kwdef struct FitParamSpec
    name::Symbol
    lo::Float64
    hi::Float64
    initial::Float64
end

Base.@kwdef struct FitVectorBundle
    field_T::Vector{Float64} = Float64[]
    qtag::Vector{String} = String[]
    energy_meV::Vector{Float64} = Float64[]
    data_intensity::Vector{Float64} = Float64[]
    data_error::Vector{Float64} = Float64[]
    model_dispersive::Vector{Float64} = Float64[]
    model_nondispersive::Vector{Float64} = Float64[]
    weight::Vector{Float64} = Float64[]
end

function _safe_initial(x::Real, lo::Real, hi::Real)
    width = Float64(hi - lo)
    width > 0 || error("Bad parameter bounds: lo=$lo hi=$hi")
    epsx = max(1e-9 * width, 1e-12)
    return clamp(Float64(x), Float64(lo) + epsx, Float64(hi) - epsx)
end

function _invlogit_stable(u::Real)
    x = Float64(u)
    if x >= 0
        z = exp(-x)
        return 1.0 / (1.0 + z)
    else
        z = exp(x)
        return z / (1.0 + z)
    end
end

function _logit(x::Real)
    y = clamp(Float64(x), 1e-12, 1.0 - 1e-12)
    return log(y / (1.0 - y))
end

function unconstrained_initial(specs::Vector{FitParamSpec})
    u0 = Float64[]
    for s in specs
        x0 = _safe_initial(s.initial, s.lo, s.hi)
        frac = (x0 - s.lo) / (s.hi - s.lo)
        push!(u0, _logit(frac))
    end
    return u0
end

function unpack_unconstrained(u::AbstractVector{<:Real}, specs::Vector{FitParamSpec})
    length(u) == length(specs) || error("Parameter vector length mismatch")
    p = Dict{Symbol, Float64}()
    for (ui, s) in zip(u, specs)
        p[s.name] = s.lo + (s.hi - s.lo) * _invlogit_stable(ui)
    end
    return p
end

function default_fit_param_specs(; scale_mode::Symbol=:analytic,
                                  initial_scale::Real=1.0,
                                  initial_second_kernel_relative_intensity::Real=1.0,
                                  initial_gzz2::Real=2.88,
                                  initial_sigma_gzz2::Real=0.25,
                                  initial_gperp::Union{Nothing,Real}=nothing,
                                  initial_gperp2::Union{Nothing,Real}=nothing,
                                  gzz_bounds=(1.0, 8.0),
                                  J1_bounds=(0.0, 0.60),
                                  J2_bounds=(0.0, 0.40),
                                  sigma_gzz_bounds=(0.0, 1.5),
                                  sigma_J1_bounds=(0.0, 2.0),
                                  sigma_J2_bounds=(0.0, 2.0),
                                  gzz2_bounds=(0.5, 8.0),
                                  sigma_gzz2_bounds=(0.0, 2.0),
                                  gperp_bounds=(0.0, 8.0),
                                  gperp2_bounds=(0.0, 8.0),
                                  log10_second_kernel_relative_intensity_bounds=(-4.0, 4.0),
                                  log10_scale_bounds=(-12.0, 12.0))
    pars = demo_defaults()
    specs = FitParamSpec[
        FitParamSpec(:gzz,       Float64(gzz_bounds[1]),       Float64(gzz_bounds[2]),       pars.model.gzz),
        FitParamSpec(:J1_meV,    Float64(J1_bounds[1]),        Float64(J1_bounds[2]),        pars.model.J1_meV),
        FitParamSpec(:J2_meV,    Float64(J2_bounds[1]),        Float64(J2_bounds[2]),        pars.model.J2_meV),
        FitParamSpec(:sigma_gzz, Float64(sigma_gzz_bounds[1]), Float64(sigma_gzz_bounds[2]), pars.disorder.sigma_gzz),
        FitParamSpec(:sigma_J1,  Float64(sigma_J1_bounds[1]),  Float64(sigma_J1_bounds[2]),  pars.disorder.sigma_J1),
        FitParamSpec(:sigma_J2,  Float64(sigma_J2_bounds[1]),  Float64(sigma_J2_bounds[2]),  pars.disorder.sigma_J2),
        FitParamSpec(:gzz2,      Float64(gzz2_bounds[1]),      Float64(gzz2_bounds[2]),      Float64(initial_gzz2)),
        FitParamSpec(:sigma_gzz2,Float64(sigma_gzz2_bounds[1]),Float64(sigma_gzz2_bounds[2]),Float64(initial_sigma_gzz2)),
        FitParamSpec(:gperp,     Float64(gperp_bounds[1]),     Float64(gperp_bounds[2]),     Float64(initial_gperp === nothing ? pars.model.gzz : initial_gperp)),
        FitParamSpec(:gperp2,    Float64(gperp2_bounds[1]),    Float64(gperp2_bounds[2]),    Float64(initial_gperp2 === nothing ? initial_gzz2 : initial_gperp2)),
        FitParamSpec(:log10_second_kernel_relative_intensity,
                     Float64(log10_second_kernel_relative_intensity_bounds[1]),
                     Float64(log10_second_kernel_relative_intensity_bounds[2]),
                     log10(Float64(initial_second_kernel_relative_intensity))),
    ]

    if scale_mode == :parameter
        push!(specs, FitParamSpec(:log10_scale,
                                  Float64(log10_scale_bounds[1]),
                                  Float64(log10_scale_bounds[2]),
                                  log10(Float64(initial_scale))))
    elseif scale_mode != :analytic
        error("scale_mode must be :analytic or :parameter")
    end

    return specs
end

function print_fit_param_specs(specs::Vector{FitParamSpec})
    println("Fit parameter bounds")
    println("--------------------")
    for s in specs
        println(@sprintf("%-42s initial=% .8g   bounds=[% .8g, % .8g]",
                         String(s.name), s.initial, s.lo, s.hi))
    end
end

function second_kernel_relative_intensity(p::Dict{Symbol,Float64})
    return 10.0 ^ p[:log10_second_kernel_relative_intensity]
end

# Neutron intensity prefactors for an effective S=1/2 Yb doublet.
# gzz controls Zeeman splitting and field-parallel magnetization, while the
# one-magnon / single-ion INS intensity is proportional to transverse matrix
# elements. In this phenomenological effective-spin version, the transverse
# matrix element is represented by gperp^2. If these keys are absent, recover
# the older script behavior with unit intensity prefactors.
function neutron_gperp_intensity_scales(p::Dict{Symbol,Float64})
    gperp_disp = haskey(p, :gperp)  ? max(p[:gperp],  0.0) : 1.0
    gperp_flat = haskey(p, :gperp2) ? max(p[:gperp2], 0.0) : 1.0
    return (; dispersive = gperp_disp^2, nondispersive = gperp_flat^2)
end

function print_fit_params(p::Dict{Symbol,Float64}; scale::Union{Nothing,Float64}=nothing)
    println("Fit parameters")
    println("--------------")
    for nm in [:gzz, :J1_meV, :J2_meV, :sigma_gzz, :sigma_J1, :sigma_J2,
               :gzz2, :sigma_gzz2, :gperp, :gperp2, :log10_second_kernel_relative_intensity,
               :log10_scale]
        if haskey(p, nm)
            println(@sprintf("%-42s % .10g", String(nm), p[nm]))
        end
    end
    if haskey(p, :log10_second_kernel_relative_intensity)
        println(@sprintf("%-42s % .10g", "second_kernel_relative_intensity", second_kernel_relative_intensity(p)))
    end
    if scale !== nothing
        println(@sprintf("%-42s % .10g", "global_scale", scale))
        println(@sprintf("%-42s % .10g", "second_kernel_effective_scale", Float64(scale) * second_kernel_relative_intensity(p)))
    end
end

function _dispersive_model_params_from_fit_dict(p::Dict{Symbol,Float64}, B_T::Real; S::Real=0.5)
    return ModelParams(
        B_T = Float64(B_T),
        gzz = p[:gzz],
        J1_meV = p[:J1_meV],
        J2_meV = p[:J2_meV],
        S = Float64(S),
    )
end

function _dispersive_disorder_params_from_fit_dict(p::Dict{Symbol,Float64};
                                                   J_units::Symbol=:fractional,
                                                   correlate_J1_J2::Bool=false)
    return DisorderParams(
        enabled = true,
        sigma_gzz = max(p[:sigma_gzz], 0.0),
        sigma_J1 = max(p[:sigma_J1], 0.0),
        sigma_J2 = max(p[:sigma_J2], 0.0),
        J_units = J_units,
        correlate_J1_J2 = correlate_J1_J2,
    )
end

function _nondispersive_model_params_from_fit_dict(p::Dict{Symbol,Float64}, B_T::Real; S::Real=0.5)
    # Effective S=1/2 kernel with no exchange terms. The energy is simply
    # E = gzz2 * μB * B. S does not affect the energy when J1=J2=0, but keeping
    # S here makes the kernel explicit and parallel to the dispersive model.
    return ModelParams(
        B_T = Float64(B_T),
        gzz = p[:gzz2],
        J1_meV = 0.0,
        J2_meV = 0.0,
        S = Float64(S),
    )
end

function _nondispersive_disorder_params_from_fit_dict(p::Dict{Symbol,Float64})
    return DisorderParams(
        enabled = true,
        sigma_gzz = max(p[:sigma_gzz2], 0.0),
        sigma_J1 = 0.0,
        sigma_J2 = 0.0,
        J_units = :absolute_meV,
        correlate_J1_J2 = false,
    )
end

function simulate_two_kernel_model_fields(p::Dict{Symbol,Float64};
                                          fields_T::AbstractVector{<:Real}=[9.0],
                                          cuts::Vector{CutSpec1D}=default_cuts_1d(),
                                          lattice::LatticeParams=demo_defaults().lattice,
                                          resolution::ResolutionParams=demo_defaults().resolution,
                                          n_samples_per_cut::Int=100_000,
                                          seed::Int=2026,
                                          S::Real=0.5,
                                          J_units::Symbol=:fractional,
                                          correlate_J1_J2::Bool=false,
                                          use_form_factor::Bool=true,
                                          include_j2_formfactor::Bool=true,
                                          include_kfki::Bool=true,
                                          polarization::Symbol=:transverse_c)
    disorder_disp = _dispersive_disorder_params_from_fit_dict(p;
        J_units=J_units,
        correlate_J1_J2=correlate_J1_J2,
    )
    disorder_flat = _nondispersive_disorder_params_from_fit_dict(p)

    out = Dict{Float64, Any}()
    for B in fields_T
        Bf = Float64(B)
        model_disp = _dispersive_model_params_from_fit_dict(p, Bf; S=S)
        model_flat = _nondispersive_model_params_from_fit_dict(p, Bf; S=S)

        # Fixed common random numbers per field make the objective deterministic.
        # The two kernels use independent but reproducible streams.
        field_seed = seed + round(Int, 1000 * Bf)
        rng_disp = MersenneTwister(field_seed)
        rng_flat = MersenneTwister(field_seed + 7919)

        gperp_scales = neutron_gperp_intensity_scales(p)

        dispersive = simulate_all_cuts_1d_stratified(
            model_disp;
            cuts=cuts,
            lattice=lattice,
            disorder=disorder_disp,
            resolution=resolution,
            rng=rng_disp,
            n_samples_per_cut=n_samples_per_cut,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
            intensity_scale=gperp_scales.dispersive,
        )

        nondispersive = simulate_all_cuts_1d_stratified(
            model_flat;
            cuts=cuts,
            lattice=lattice,
            disorder=disorder_flat,
            resolution=resolution,
            rng=rng_flat,
            n_samples_per_cut=n_samples_per_cut,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
            intensity_scale=gperp_scales.nondispersive,
        )

        out[Bf] = (dispersive=dispersive, nondispersive=nondispersive)
    end
    return out
end

function _fit_window_mask(E::AbstractVector{<:Real}, qtag::String, fit_windows_by_q)
    windows = fit_windows_by_q === nothing ? nothing : get(fit_windows_by_q, qtag, nothing)
    return _window_mask(E, windows)
end

function _combined_model_vector(bundle::FitVectorBundle, r2::Real)
    return bundle.model_dispersive .+ Float64(r2) .* bundle.model_nondispersive
end

function build_fit_vector_bundle(data_scans,
                                 model_results_by_field;
                                 fields_T::AbstractVector{<:Real}=[9.0],
                                 qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                                 qtag_to_modelkey=OVERLAY_QTAG_TO_MODELKEY,
                                 fit_windows_by_q=YZGO_FIT_WINDOWS_1D,
                                 use_errors::Bool=true,
                                 error_floor::Real=0.0)
    bundle = FitVectorBundle()

    for B in fields_T
        Bf = Float64(B)
        if !haskey(model_results_by_field, Bf)
            @warn "Skipping field $Bf T because model results are absent."
            continue
        end
        model_components = model_results_by_field[Bf]
        model_disp_results = model_components.dispersive
        model_flat_results = model_components.nondispersive

        for qtag in qtags
            if !haskey(data_scans, qtag)
                @warn "Skipping qtag $qtag because it is absent from data_scans."
                continue
            end
            if !haskey(data_scans[qtag], Bf)
                @warn "Skipping qtag $qtag at $Bf T because that field is absent from data_scans."
                continue
            end

            scan = data_scans[qtag][Bf]
            mres_disp = _model_result_for_qtag(model_disp_results, qtag; qtag_to_modelkey=qtag_to_modelkey)
            mres_flat = _model_result_for_qtag(model_flat_results, qtag; qtag_to_modelkey=qtag_to_modelkey)
            model_disp_on_data = interp1_linear_local(mres_disp.E_centers_meV, mres_disp.intensity, scan.energy)
            model_flat_on_data = interp1_linear_local(mres_flat.E_centers_meV, mres_flat.intensity, scan.energy)

            fit_mask = isfinite.(scan.energy) .&
                       isfinite.(scan.intensity) .&
                       isfinite.(model_disp_on_data) .&
                       isfinite.(model_flat_on_data) .&
                       _fit_window_mask(scan.energy, qtag, fit_windows_by_q)

            if use_errors
                fit_mask .&= isfinite.(scan.error) .& (scan.error .> 0)
            end

            nuse = count(fit_mask)
            @info "Fit mask points" field_T=Bf qtag=qtag windows=(fit_windows_by_q === nothing ? nothing : get(fit_windows_by_q, qtag, nothing)) nuse=nuse
            if nuse < 3
                @warn "Very few usable fit points for qtag=$qtag, field=$(Bf) T; count=$nuse."
            end

            for i in findall(fit_mask)
                err = use_errors ? sqrt(Float64(scan.error[i])^2 + Float64(error_floor)^2) : 1.0
                err = max(err, eps(Float64))
                push!(bundle.field_T, Bf)
                push!(bundle.qtag, qtag)
                push!(bundle.energy_meV, Float64(scan.energy[i]))
                push!(bundle.data_intensity, Float64(scan.intensity[i]))
                push!(bundle.data_error, err)
                push!(bundle.model_dispersive, Float64(model_disp_on_data[i]))
                push!(bundle.model_nondispersive, Float64(model_flat_on_data[i]))
                push!(bundle.weight, 1.0 / err^2)
            end
        end
    end

    return bundle
end

function best_global_scale(bundle::FitVectorBundle, r2::Real; positive_only::Bool=true)
    y = bundle.data_intensity
    m = _combined_model_vector(bundle, r2)
    w = bundle.weight
    length(y) >= 3 || return NaN
    denom = sum(w .* m .* m)
    denom > 0 || return NaN
    s = sum(w .* y .* m) / denom
    return positive_only ? max(s, 0.0) : s
end

function chisq_for_bundle(bundle::FitVectorBundle, scale::Real, r2::Real; nfree::Int=1, reduced::Bool=true)
    y = bundle.data_intensity
    m = _combined_model_vector(bundle, r2)
    w = bundle.weight
    length(y) >= 3 || return Inf
    s = Float64(scale)
    chi2 = sum(w .* (y .- s .* m).^2)
    if reduced
        dof = max(length(y) - nfree, 1)
        return chi2 / dof
    else
        return chi2
    end
end

function objective_for_unconstrained(u::AbstractVector{<:Real},
                                     specs::Vector{FitParamSpec},
                                     data_scans;
                                     fields_T::AbstractVector{<:Real}=[9.0],
                                     qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                                     fit_windows_by_q=YZGO_FIT_WINDOWS_1D,
                                     scale_mode::Symbol=:analytic,
                                     n_samples_per_cut::Int=100_000,
                                     seed::Int=2026,
                                     use_errors::Bool=true,
                                     error_floor::Real=0.0,
                                     lattice::LatticeParams=demo_defaults().lattice,
                                     resolution::ResolutionParams=demo_defaults().resolution,
                                     cuts::Vector{CutSpec1D}=default_cuts_1d(),
                                     S::Real=0.5,
                                     J_units::Symbol=:fractional,
                                     correlate_J1_J2::Bool=false,
                                     use_form_factor::Bool=true,
                                     include_j2_formfactor::Bool=true,
                                     include_kfki::Bool=true,
                                     polarization::Symbol=:transverse_c,
                                     penalty::Real=1e30)
    try
        p = unpack_unconstrained(u, specs)
        r2 = second_kernel_relative_intensity(p)

        model_results_by_field = simulate_two_kernel_model_fields(p;
            fields_T=fields_T,
            cuts=cuts,
            lattice=lattice,
            resolution=resolution,
            n_samples_per_cut=n_samples_per_cut,
            seed=seed,
            S=S,
            J_units=J_units,
            correlate_J1_J2=correlate_J1_J2,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
        )

        bundle = build_fit_vector_bundle(data_scans, model_results_by_field;
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            use_errors=use_errors,
            error_floor=error_floor,
        )

        if length(bundle.data_intensity) < 6
            return penalty
        end

        scale = if scale_mode == :analytic
            best_global_scale(bundle, r2)
        elseif scale_mode == :parameter
            10.0 ^ p[:log10_scale]
        else
            error("scale_mode must be :analytic or :parameter")
        end

        if !isfinite(scale) || scale < 0
            return penalty
        end

        nfree = length(specs) + (scale_mode == :analytic ? 1 : 0)
        val = chisq_for_bundle(bundle, scale, r2; nfree=nfree, reduced=true)
        return isfinite(val) ? val : penalty
    catch err
        @warn "Objective evaluation failed" exception=(err, catch_backtrace())
        return penalty
    end
end

function write_fit_points_csv(path::AbstractString, bundle::FitVectorBundle, scale::Real, r2::Real)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "field_T,qtag,energy_meV,data_intensity,data_error,model_dispersive_unscaled,model_nondispersive_unscaled,second_kernel_relative_intensity,model_combined_unscaled,global_scale,model_scaled,residual,weight")
        for i in eachindex(bundle.energy_meV)
            model_combined = bundle.model_dispersive[i] + Float64(r2) * bundle.model_nondispersive[i]
            model_scaled = Float64(scale) * model_combined
            residual = bundle.data_intensity[i] - model_scaled
            println(io, join((
                bundle.field_T[i],
                bundle.qtag[i],
                bundle.energy_meV[i],
                bundle.data_intensity[i],
                bundle.data_error[i],
                bundle.model_dispersive[i],
                bundle.model_nondispersive[i],
                Float64(r2),
                model_combined,
                Float64(scale),
                model_scaled,
                residual,
                bundle.weight[i],
            ), ","))
        end
    end
    return path
end

function write_scaled_model_csv(path::AbstractString,
                                model_results_by_field;
                                scale::Real,
                                r2::Real,
                                qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                                qtag_to_modelkey=OVERLAY_QTAG_TO_MODELKEY)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "field_T,qtag,model_key,energy_meV,model_dispersive_unscaled,model_nondispersive_unscaled,second_kernel_relative_intensity,model_combined_unscaled,global_scale,model_dispersive_scaled,model_nondispersive_scaled,model_total_scaled")
        for B in sort(collect(keys(model_results_by_field)))
            model_components = model_results_by_field[B]
            model_disp_results = model_components.dispersive
            model_flat_results = model_components.nondispersive
            for qtag in qtags
                model_key = qtag_to_modelkey[qtag]
                if !haskey(model_disp_results, model_key) || !haskey(model_flat_results, model_key)
                    @warn "No model result for qtag=$qtag key=$model_key"
                    continue
                end
                mres_disp = model_disp_results[model_key]
                mres_flat = model_flat_results[model_key]
                for i in eachindex(mres_disp.E_centers_meV)
                    Ed = mres_disp.E_centers_meV[i]
                    Ifloat = mres_disp.intensity[i]
                    Iflat = mres_flat.intensity[i]
                    Icomb = Ifloat + Float64(r2) * Iflat
                    println(io, join((
                        B,
                        qtag,
                        model_key,
                        Ed,
                        Ifloat,
                        Iflat,
                        Float64(r2),
                        Icomb,
                        Float64(scale),
                        Float64(scale) * Ifloat,
                        Float64(scale) * Float64(r2) * Iflat,
                        Float64(scale) * Icomb,
                    ), ","))
                end
            end
        end
    end
    return path
end

function write_fit_summary(path::AbstractString;
                           specs::Vector{FitParamSpec},
                           params::Dict{Symbol,Float64},
                           scale::Real,
                           r2::Real,
                           chi2_reduced::Real,
                           npoints::Int,
                           nfree::Int,
                           fields_T,
                           qtags,
                           fit_windows_by_q,
                           n_samples_per_cut::Int,
                           final_n_samples_per_cut::Int,
                           seed::Int,
                           optimizer_result=nothing)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "YZGO 1D experimental-data fit summary: two-kernel model")
        println(io, "Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io)
        println(io, "Model:")
        println(io, "  model_scaled = global_scale * (dispersive + second_kernel_relative_intensity * nondispersive)")
        println(io, "  nondispersive kernel has J1=J2=0 and fitted gzz2, sigma_gzz2")
        println(io, "  no qtag-specific relative scale is used")
        println(io)
        println(io, "Fields fitted (T): ", collect(fields_T))
        println(io, "Q cuts fitted: ", collect(qtags))
        println(io, "Fit windows by qtag:")
        for qtag in qtags
            println(io, "  ", qtag, " => ", fit_windows_by_q === nothing ? nothing : get(fit_windows_by_q, qtag, nothing))
        end
        println(io)
        println(io, "Monte Carlo samples per cut during optimization: ", n_samples_per_cut)
        println(io, "Monte Carlo samples per cut for final model: ", final_n_samples_per_cut)
        println(io, "Common-random-number seed: ", seed)
        println(io)
        println(io, "Parameter bounds and initials:")
        for s in specs
            println(io, @sprintf("  %-42s initial=% .10g bounds=[% .10g, % .10g]",
                                 String(s.name), s.initial, s.lo, s.hi))
        end
        println(io)
        println(io, "Best-fit parameters:")
        for nm in [:gzz, :J1_meV, :J2_meV, :sigma_gzz, :sigma_J1, :sigma_J2,
                   :gzz2, :sigma_gzz2, :gperp, :gperp2, :log10_second_kernel_relative_intensity,
                   :log10_scale]
            if haskey(params, nm)
                println(io, @sprintf("  %-42s % .12g", String(nm), params[nm]))
            end
        end
        println(io, @sprintf("  %-42s % .12g", "second_kernel_relative_intensity", r2))
        println(io, @sprintf("  %-42s % .12g", "global_scale", scale))
        println(io, @sprintf("  %-42s % .12g", "second_kernel_effective_scale", Float64(scale) * Float64(r2)))
        println(io)
        println(io, @sprintf("Reduced chi^2: %.12g", chi2_reduced))
        println(io, "N fit points: ", npoints)
        println(io, "N free parameters counted in reduced chi^2: ", nfree)
        if optimizer_result !== nothing
            println(io)
            println(io, "Optimizer summary:")
            println(io, optimizer_result)
        end
    end
    return path
end

function plot_yzgo_fit_overlay_1d(data_scans,
                                  model_results_by_field;
                                  scale::Real,
                                  r2::Real,
                                  fields_T::AbstractVector{<:Real}=[9.0],
                                  qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                                  qtag_to_modelkey=OVERLAY_QTAG_TO_MODELKEY,
                                  fit_windows_by_q=YZGO_FIT_WINDOWS_1D,
                                  ylims::Union{Nothing,Tuple{Float64,Float64}}=YZGO_FIT_PLOT_YLIMS,
                                  figure_title::AbstractString="YZGO 1D two-kernel fit: data vs model",
                                  outpath::Union{Nothing,AbstractString}=nothing,
                                  save_png::Bool=true,
                                  display_fig::Bool=false)
    fig = Figure(size=(1450, 360 * length(fields_T) + 100))
    Label(fig[0, :], figure_title; fontsize=20)

    for (irow, B) in enumerate(fields_T)
        Bf = Float64(B)
        haskey(model_results_by_field, Bf) || continue
        model_components = model_results_by_field[Bf]
        model_disp_results = model_components.dispersive
        model_flat_results = model_components.nondispersive

        for (icol, qtag) in enumerate(qtags)
            if !(haskey(data_scans, qtag) && haskey(data_scans[qtag], Bf))
                continue
            end

            scan = data_scans[qtag][Bf]
            title = get(OVERLAY_QTAG_TITLES, qtag, qtag)
            ax = Axis(fig[irow, icol];
                title=@sprintf("%.0f T, %s", Bf, title),
                xlabel="Delta E (meV)",
                ylabel=icol == 1 ? "Intensity - background" : "",
            )

            scatter!(ax, scan.energy, scan.intensity; label="data")
            errorbars!(ax, scan.energy, scan.intensity, scan.error; whiskerwidth=4)

            mres_disp = _model_result_for_qtag(model_disp_results, qtag; qtag_to_modelkey=qtag_to_modelkey)
            mres_flat = _model_result_for_qtag(model_flat_results, qtag; qtag_to_modelkey=qtag_to_modelkey)
            E = mres_disp.E_centers_meV
            y_disp = Float64(scale) .* mres_disp.intensity
            y_flat = Float64(scale) .* Float64(r2) .* mres_flat.intensity
            y_total = y_disp .+ y_flat

            lines!(ax, E, y_disp; label="dispersive", linewidth=2, linestyle=:dash)
            lines!(ax, E, y_flat; label="non-dispersive", linewidth=2, linestyle=:dot)
            lines!(ax, E, y_total; label="total", linewidth=3)

            windows = fit_windows_by_q === nothing ? Tuple{Float64,Float64}[] : get(fit_windows_by_q, qtag, Tuple{Float64,Float64}[])
            for (lo, hi) in windows
                vlines!(ax, [lo, hi]; linestyle=:dash, linewidth=1)
            end

            if hasproperty(scan, :meta)
                xlims!(ax, scan.meta.de_range...)
            end
            if ylims !== nothing
                ylims!(ax, ylims...)
            end
            axislegend(ax; position=:rt, framevisible=false)
        end
    end

    if save_png && outpath !== nothing
        mkpath(dirname(outpath))
        save(outpath, fig)
        println("Saved fit overlay figure to: ", outpath)
    end
    if display_fig
        display(fig)
    end
    return fig
end

function run_yzgo_fit_1d(; base_dir::AbstractString=BASEDIR,
                         outdir::AbstractString=joinpath(base_dir, "fit_model_1d_two_kernel"),
                         fields_T::AbstractVector{<:Real}=[9.0],
                         qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                         fit_windows_by_q=YZGO_FIT_WINDOWS_1D,
                         data_mode::Symbol=:tail_bgsub,
                         scale_mode::Symbol=:analytic,
                         plot_ylims::Union{Nothing,Tuple{Float64,Float64}}=YZGO_FIT_PLOT_YLIMS,
                         specs::Union{Nothing,Vector{FitParamSpec}}=nothing,
                         maxiters::Int=80,
                         n_samples_per_cut::Int=100_000,
                         final_n_samples_per_cut::Int=500_000,
                         seed::Int=2026,
                         use_errors::Bool=true,
                         error_floor::Real=0.0,
                         show_trace::Bool=true,
                         trace_every::Int=5,
                         make_plots::Bool=true,
                         display_figures::Bool=true,
                         lattice::LatticeParams=demo_defaults().lattice,
                         resolution::ResolutionParams=demo_defaults().resolution,
                         cuts::Vector{CutSpec1D}=default_cuts_1d(),
                         S::Real=0.5,
                         J_units::Symbol=:fractional,
                         correlate_J1_J2::Bool=false,
                         use_form_factor::Bool=true,
                         include_j2_formfactor::Bool=true,
                         include_kfki::Bool=true,
                         polarization::Symbol=:transverse_c,
                         background_kwargs...)
    mkpath(outdir)

    specs2 = specs === nothing ? default_fit_param_specs(scale_mode=scale_mode) : specs
    print_fit_param_specs(specs2)

    println()
    println("Loading experimental scans from: ", base_dir)
    data_scans, background_models, scans_raw = load_fit_data_1d(;
        base_dir=base_dir,
        data_mode=data_mode,
        background_kwargs...
    )

    println("Using data mode: ", data_mode)
    println("Fitting fields: ", collect(fields_T))
    println("Fitting q cuts: ", collect(qtags))
    println("Fit windows: ", fit_windows_by_q)
    println("Second kernel: non-dispersive S=1/2 with J1=J2=0")
    println("Second-kernel initial guess: gzz2=2.88, sigma_gzz2=0.25")
    println("No qtag-specific relative scale is used in this script.")
    println("Plot y-limits: ", plot_ylims)
    println("Optimization samples per cut: ", n_samples_per_cut)
    println("Final samples per cut: ", final_n_samples_per_cut)
    println()

    u0 = unconstrained_initial(specs2)
    eval_counter = Ref(0)

    obj = function(u)
        eval_counter[] += 1
        val = objective_for_unconstrained(u, specs2, data_scans;
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            scale_mode=scale_mode,
            n_samples_per_cut=n_samples_per_cut,
            seed=seed,
            use_errors=use_errors,
            error_floor=error_floor,
            lattice=lattice,
            resolution=resolution,
            cuts=cuts,
            S=S,
            J_units=J_units,
            correlate_J1_J2=correlate_J1_J2,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
        )

        if trace_every > 0 && (eval_counter[] == 1 || eval_counter[] % trace_every == 0)
            pnow = unpack_unconstrained(u, specs2)
            r2now = second_kernel_relative_intensity(pnow)
            @printf("eval %4d  redchi2 = %.8g   gzz=%.5g  J1=%.5g  J2=%.5g  sigg=%.5g  sigJ1=%.5g  sigJ2=%.5g  gzz2=%.5g  sigg2=%.5g  r2=%.5g\n",
                    eval_counter[], val,
                    pnow[:gzz], pnow[:J1_meV], pnow[:J2_meV],
                    pnow[:sigma_gzz], pnow[:sigma_J1], pnow[:sigma_J2],
                    pnow[:gzz2], pnow[:sigma_gzz2], r2now)
        end

        return val
    end

    println("Starting Optim.jl Nelder-Mead fit in bounded-transformed parameters...")
    optres = optimize(obj, u0, NelderMead(), Optim.Options(iterations=maxiters, show_trace=show_trace))

    ubest = Optim.minimizer(optres)
    pbest = unpack_unconstrained(ubest, specs2)
    r2_final = second_kernel_relative_intensity(pbest)

    println()
    println("Best optimizer result:")
    println(optres)
    println()

    println("Rerunning final two-kernel model with final_n_samples_per_cut = ", final_n_samples_per_cut)
    model_final = simulate_two_kernel_model_fields(pbest;
        fields_T=fields_T,
        cuts=cuts,
        lattice=lattice,
        resolution=resolution,
        n_samples_per_cut=final_n_samples_per_cut,
        seed=seed,
        S=S,
        J_units=J_units,
        correlate_J1_J2=correlate_J1_J2,
        use_form_factor=use_form_factor,
        include_j2_formfactor=include_j2_formfactor,
        include_kfki=include_kfki,
        polarization=polarization,
    )

    bundle_final = build_fit_vector_bundle(data_scans, model_final;
        fields_T=fields_T,
        qtags=qtags,
        fit_windows_by_q=fit_windows_by_q,
        use_errors=use_errors,
        error_floor=error_floor,
    )

    scale_final = if scale_mode == :analytic
        best_global_scale(bundle_final, r2_final)
    elseif scale_mode == :parameter
        10.0 ^ pbest[:log10_scale]
    else
        error("scale_mode must be :analytic or :parameter")
    end

    nfree = length(specs2) + (scale_mode == :analytic ? 1 : 0)
    redchi2_final = chisq_for_bundle(bundle_final, scale_final, r2_final; nfree=nfree, reduced=true)

    println()
    print_fit_params(pbest; scale=scale_final)
    println(@sprintf("Reduced chi^2 = %.8g from %d fit points", redchi2_final, length(bundle_final.data_intensity)))
    println()

    summary_path = joinpath(outdir, "YZGO_1d_two_kernel_fit_summary.txt")
    fit_points_path = joinpath(outdir, "YZGO_1d_two_kernel_fit_points.csv")
    model_csv_path = joinpath(outdir, "YZGO_1d_two_kernel_scaled_model.csv")
    overlay_path = joinpath(outdir, "YZGO_1d_two_kernel_fit_overlay.png")

    write_fit_summary(summary_path;
        specs=specs2,
        params=pbest,
        scale=scale_final,
        r2=r2_final,
        chi2_reduced=redchi2_final,
        npoints=length(bundle_final.data_intensity),
        nfree=nfree,
        fields_T=fields_T,
        qtags=qtags,
        fit_windows_by_q=fit_windows_by_q,
        n_samples_per_cut=n_samples_per_cut,
        final_n_samples_per_cut=final_n_samples_per_cut,
        seed=seed,
        optimizer_result=optres,
    )
    write_fit_points_csv(fit_points_path, bundle_final, scale_final, r2_final)
    write_scaled_model_csv(model_csv_path, model_final;
        scale=scale_final,
        r2=r2_final,
        qtags=qtags,
    )

    fig = nothing
    if make_plots
        fig = plot_yzgo_fit_overlay_1d(data_scans, model_final;
            scale=scale_final,
            r2=r2_final,
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            ylims=plot_ylims,
            outpath=overlay_path,
            save_png=true,
            display_fig=display_figures,
        )
    end

    println("Wrote:")
    println("  ", summary_path)
    println("  ", fit_points_path)
    println("  ", model_csv_path)
    if make_plots
        println("  ", overlay_path)
    end

    return (;
        optimizer_result=optres,
        params=pbest,
        scale=scale_final,
        second_kernel_relative_intensity=r2_final,
        second_kernel_effective_scale=scale_final * r2_final,
        reduced_chisq=redchi2_final,
        nfit_points=length(bundle_final.data_intensity),
        specs=specs2,
        data_scans=data_scans,
        scans_raw=scans_raw,
        background_models=background_models,
        model_results=model_final,
        fit_vectors=bundle_final,
        figure=fig,
        paths=(summary=summary_path, fit_points=fit_points_path, model=model_csv_path, overlay=overlay_path),
    )
end


end # module YZGONeutronFit

module YZGOMagnetizationFit
# YZGO_fit_full_field_magnetization_impurity.jl
#
# Fit the full-field digitized magnetization curve of YbZn2GaO5 / YZGO to a
# two-component onsite-polarization model:
#
#   1. Dispersive/background Yb component with J1/J2 disorder and gzz disorder.
#      This is the same fast analytical component used in the field-polarized
#      model.  sigma_J2 is fixed equal to sigma_J1.
#
#   2. Non-dispersive S = 1/2 impurity-like component with no exchange terms,
#      but with a mean gzz2 and Gaussian gzz2 disorder.
#
# The two components have independent fitted scale factors:
#
#   M_total(B) = scale_disp * (M_disp(B) + chi_vv_muB_per_T * B)
#              + scale_imp  *  M_imp(B)
#
# where chi_vv_muB_per_T is a phenomenological Van Vleck slope in units of
# mu_B / T / Yb before multiplication by scale_disp.
#
# Usage from Julia:
#     include("YZGO_fit_full_field_magnetization_impurity.jl")
#     out = fit_full_field_magnetization_impurity()
#
# Or run directly:
#     julia YZGO_fit_full_field_magnetization_impurity.jl

using Random
using Statistics
using Printf

const MU_B_MEV_PER_T = 5.7883818060e-2
const KB_MEV_PER_K   = 8.617333262e-2

const HAS_GLMAKIE = let ok = false
    try
        @eval using GLMakie
        ok = true
    catch
        ok = false
    end
    ok
end

Base.@kwdef struct FitParams
    # Independent scale factors for the two modeled constituents.
    scale_disp::Float64 = 0.5
    scale_imp::Float64 = 0.05

    # Phenomenological Van Vleck term, attached to the dispersive/Yb component.
    chi_vv_muB_per_T::Float64 = 0.07

    # Dispersive component disorder and mean parameters.
    sigma_gzz::Float64 = 1.0 / 3.0
    sigma_J1::Float64 = 1.0 / 3.0
    gzz::Float64 = 3.5
    J1_meV::Float64 = 0.177824
    J2_meV::Float64 = 0.00973612

    # Non-dispersive impurity-like S = 1/2 component.
    gzz2::Float64 = 2.88
    sigma_gzz2::Float64 = 0.5
end

function vector_from_params(p::FitParams)
    return [
        p.scale_disp,
        p.scale_imp,
        p.chi_vv_muB_per_T,
        p.sigma_gzz,
        p.sigma_J1,
        p.gzz,
        p.J1_meV,
        p.J2_meV,
        p.gzz2,
        p.sigma_gzz2,
    ]
end

function params_from_vector(x::AbstractVector{<:Real})
    return FitParams(
        scale_disp        = Float64(x[1]),
        scale_imp         = Float64(x[2]),
        chi_vv_muB_per_T  = Float64(x[3]),
        sigma_gzz         = Float64(x[4]),
        sigma_J1          = Float64(x[5]),
        gzz               = Float64(x[6]),
        J1_meV            = Float64(x[7]),
        J2_meV            = Float64(x[8]),
        gzz2              = Float64(x[9]),
        sigma_gzz2        = Float64(x[10]),
    )
end

function clamp_vector(x, lower, upper)
    y = similar(Float64.(x))
    for i in eachindex(x)
        y[i] = clamp(Float64(x[i]), lower[i], upper[i])
    end
    return y
end

function safe_ratio(a::Real, b::Real)
    return abs(b) > eps(Float64) ? Float64(a) / Float64(b) : NaN
end

function component_fraction(p::FitParams)
    denom = p.scale_disp + p.scale_imp
    if abs(denom) <= eps(Float64)
        return (; disp = NaN, imp = NaN)
    end
    return (; disp = p.scale_disp / denom, imp = p.scale_imp / denom)
end

function print_params(p::FitParams; prefix::AbstractString = "")
    frac = component_fraction(p)
    println(prefix * @sprintf("scale_disp       = %.8g", p.scale_disp))
    println(prefix * @sprintf("scale_imp        = %.8g", p.scale_imp))
    println(prefix * @sprintf("relative disp    = %.8g", frac.disp))
    println(prefix * @sprintf("relative imp     = %.8g", frac.imp))
    println(prefix * @sprintf("chi_vv_muB_per_T = %.8g", p.chi_vv_muB_per_T))
    println(prefix * @sprintf("sigma_gzz        = %.8g", p.sigma_gzz))
    println(prefix * @sprintf("sigma_J1_frac    = %.8g", p.sigma_J1))
    println(prefix * @sprintf("sigma_J2_frac    = %.8g  (fixed equal to sigma_J1_frac)", p.sigma_J1))
    println(prefix * @sprintf("gzz              = %.8g", p.gzz))
    println(prefix * @sprintf("J1_meV           = %.8g", p.J1_meV))
    println(prefix * @sprintf("J2_meV           = %.8g", p.J2_meV))
    println(prefix * @sprintf("J2/J1            = %.8g", safe_ratio(p.J2_meV, p.J1_meV)))
    println(prefix * @sprintf("gzz2             = %.8g", p.gzz2))
    println(prefix * @sprintf("sigma_gzz2       = %.8g", p.sigma_gzz2))
end

function write_params(io::IO, p::FitParams)
    frac = component_fraction(p)
    println(io, @sprintf("scale_disp       = %.8g", p.scale_disp))
    println(io, @sprintf("scale_imp        = %.8g", p.scale_imp))
    println(io, @sprintf("relative disp    = %.8g", frac.disp))
    println(io, @sprintf("relative imp     = %.8g", frac.imp))
    println(io, @sprintf("chi_vv_muB_per_T = %.8g", p.chi_vv_muB_per_T))
    println(io, @sprintf("sigma_gzz        = %.8g", p.sigma_gzz))
    println(io, @sprintf("sigma_J1_frac    = %.8g", p.sigma_J1))
    println(io, @sprintf("sigma_J2_frac    = %.8g  (fixed equal to sigma_J1_frac)", p.sigma_J1))
    println(io, @sprintf("gzz              = %.8g", p.gzz))
    println(io, @sprintf("J1_meV           = %.8g", p.J1_meV))
    println(io, @sprintf("J2_meV           = %.8g", p.J2_meV))
    println(io, @sprintf("J2/J1            = %.8g", safe_ratio(p.J2_meV, p.J1_meV)))
    println(io, @sprintf("gzz2             = %.8g", p.gzz2))
    println(io, @sprintf("sigma_gzz2       = %.8g", p.sigma_gzz2))
end

# -----------------------------
# CSV reader for two-column-or-more magnetization data
# -----------------------------

function read_experimental_magnetization_csv(filename::AbstractString)
    if !isfile(filename)
        error("Could not find experimental CSV: $(filename)")
    end

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
            if isnothing(b) || isnothing(m)
                # Skip header and non-numeric rows.
                continue
            end
            push!(B, b)
            push!(M, m)
        end
    end

    isempty(B) && error("No numeric rows found in $(filename)")

    p = sortperm(B)
    return (; B_T = B[p], M_muB_per_Yb = M[p], filename = filename)
end

function default_experiment_path()
    candidates = [
        raw"C:\Users\vdp\ORNL Dropbox\Daniel Pajerowski\YZGO\YZGO\data\YZGO_MvB_black_curve_digitized_visible.csv",
        joinpath(@__DIR__, "YZGO_MvB_black_curve_digitized_visible.csv"),
        joinpath(@__DIR__, "..", "data", "YZGO_MvB_black_curve_digitized_visible.csv"),
    ]
    for c in candidates
        if isfile(c)
            return c
        end
    end
    return candidates[1]
end

function filter_fit_window(experiment; Bmin_fit_T::Real = 0.0, Bmax_fit_T::Real = 7.0)
    B = Float64[]
    M = Float64[]
    for (b, m) in zip(experiment.B_T, experiment.M_muB_per_Yb)
        if isfinite(b) && isfinite(m) && b >= Float64(Bmin_fit_T) && b <= Float64(Bmax_fit_T)
            push!(B, b)
            push!(M, m)
        end
    end
    length(B) < 4 && error("Too few experimental points in fit window $(Bmin_fit_T) to $(Bmax_fit_T) T")
    return (; B_T = B, M_muB_per_Yb = M, filename = experiment.filename)
end

# -----------------------------
# Disorder draws and quadrature
# -----------------------------

function make_standard_normal_draws(n_samples::Int; seed::Int = 20260520)
    rng = MersenneTwister(seed)
    return (; zg = randn(rng, n_samples), z1 = randn(rng, n_samples), z2 = randn(rng, n_samples))
end

function make_normal_quadrature(n::Int = 101; zmax::Real = 5.0)
    # Lightweight deterministic quadrature for E[f(z)] with z ~ N(0,1).
    # Uses a uniformly spaced grid in z and normalized normal weights.
    n >= 5 || error("normal quadrature needs at least 5 points")
    z = collect(range(-Float64(zmax), Float64(zmax); length = n))
    w = exp.(-0.5 .* z .^ 2)
    w ./= sum(w)
    return (; z, w)
end

# -----------------------------
# Dispersive component
# -----------------------------

function exchange_Dmax_meV(J1_meV::Real, J2_meV::Real; mode::Symbol = :high_symmetry)
    J1 = Float64(J1_meV)
    J2 = Float64(J2_meV)
    if mode == :high_symmetry
        return max(0.0, 9.0 * J1, 8.0 * (J1 + J2))
    elseif mode == :K_only
        return max(0.0, 9.0 * J1)
    elseif mode == :M_only
        return max(0.0, 8.0 * (J1 + J2))
    else
        error("Unknown dmax mode $(mode). Use :high_symmetry, :K_only, or :M_only.")
    end
end

function clean_Bsat_T(p::FitParams; S::Real = 0.5, dmax_mode::Symbol = :high_symmetry)
    Dmax = exchange_Dmax_meV(p.J1_meV, p.J2_meV; mode = dmax_mode)
    return Float64(S) * Dmax / (p.gzz * MU_B_MEV_PER_T)
end

function draw_local_env_Bsat_and_msat(p::FitParams, draws;
                                      S::Real = 0.5,
                                      dmax_mode::Symbol = :high_symmetry)
    n = length(draws.zg)
    g = Vector{Float64}(undef, n)
    Bsat = Vector{Float64}(undef, n)
    msat = Vector{Float64}(undef, n)

    sigmaJ = p.sigma_J1
    sJ1_abs = sigmaJ * abs(p.J1_meV)
    sJ2_abs = sigmaJ * abs(p.J2_meV)  # sigma_J2 fixed to sigma_J1
    Sval = Float64(S)

    for i in 1:n
        gi = p.gzz + p.sigma_gzz * draws.zg[i]
        J1i = p.J1_meV + sJ1_abs * draws.z1[i]
        J2i = p.J2_meV + sJ2_abs * draws.z2[i]
        Dmaxi = exchange_Dmax_meV(J1i, J2i; mode = dmax_mode)

        g[i] = gi
        msat[i] = gi * Sval
        if isfinite(gi) && gi > 0.0 && isfinite(Dmaxi) && Dmaxi > 0.0
            Bsat[i] = Sval * Dmaxi / (gi * MU_B_MEV_PER_T)
        elseif isfinite(gi) && gi > 0.0
            Bsat[i] = 0.0
        else
            Bsat[i] = NaN
        end
    end
    return (; g, Bsat, msat)
end

function ground_magnetization_linear_saturation(B_fields_T::AbstractVector{<:Real}, env)
    # Average M_i(B) = msat_i * min(B/Bsat_i, 1).
    # Sorting makes this exact model evaluation O(N log N + NB log N).
    Bq = Float64.(B_fields_T)
    valid = isfinite.(env.Bsat) .& isfinite.(env.msat) .& (env.g .> 0.0)
    nvalid = count(valid)
    nvalid == 0 && return fill(NaN, length(Bq))

    b = env.Bsat[valid]
    m = env.msat[valid]

    zero_or_negative = b .<= 0.0
    base_sat = sum(m[zero_or_negative])

    bpos = b[.!zero_or_negative]
    mpos = m[.!zero_or_negative]

    if isempty(bpos)
        return fill(base_sat / nvalid, length(Bq))
    end

    idx = sortperm(bpos)
    bs = bpos[idx]
    ms = mpos[idx]
    prefix_m = cumsum(ms)
    slopes = ms ./ bs
    suffix_slope = reverse(cumsum(reverse(slopes)))

    out = similar(Bq)
    for (i, B) in enumerate(Bq)
        sgn = B >= 0.0 ? 1.0 : -1.0
        Ba = abs(B)
        k = searchsortedlast(bs, Ba)
        saturated_sum = k > 0 ? prefix_m[k] : 0.0
        unsat_slope_sum = k < length(bs) ? suffix_slope[k + 1] : 0.0
        out[i] = sgn * (base_sat + saturated_sum + Ba * unsat_slope_sum) / nvalid
    end
    return out
end

function ground_magnetization_smooth_crossover(B_fields_T::AbstractVector{<:Real}, env)
    # Softer onsite response:
    #     M_i(B) = msat_i * B / sqrt(B^2 + Bsat_i^2)
    # This is slower than the linear-saturation option.
    Bq = Float64.(B_fields_T)
    valid = isfinite.(env.Bsat) .& isfinite.(env.msat) .& (env.g .> 0.0)
    nvalid = count(valid)
    nvalid == 0 && return fill(NaN, length(Bq))

    b = env.Bsat[valid]
    m = env.msat[valid]
    out = zeros(Float64, length(Bq))
    for (i, B) in enumerate(Bq)
        total = 0.0
        for j in eachindex(b)
            if b[j] <= 0.0
                total += B >= 0.0 ? m[j] : -m[j]
            else
                total += m[j] * B / sqrt(B^2 + b[j]^2)
            end
        end
        out[i] = total / nvalid
    end
    return out
end

# -----------------------------
# Non-dispersive impurity component
# -----------------------------

function impurity_spinhalf_magnetization(B_fields_T::AbstractVector{<:Real}, p::FitParams, quad;
                                         temperature_K::Real = 0.4,
                                         S_imp::Real = 0.5)
    # M_imp(B) in mu_B per impurity spin, before multiplication by scale_imp.
    # For S = 1/2:
    #     M = g/2 * tanh(g mu_B B / (2 k_B T)).
    # More generally here, with S_imp defaulting to 1/2:
    #     M = g S_imp * tanh(g mu_B B S_imp / (k_B T)).
    Bq = Float64.(B_fields_T)
    T = Float64(temperature_K)
    Sval = Float64(S_imp)

    gvals_all = p.gzz2 .+ p.sigma_gzz2 .* quad.z
    valid = isfinite.(gvals_all) .& (gvals_all .> 0.0)
    if count(valid) == 0
        return fill(NaN, length(Bq))
    end
    gvals = gvals_all[valid]
    w = quad.w[valid]
    w ./= sum(w)

    out = zeros(Float64, length(Bq))
    if T <= 0.0
        for (i, B) in enumerate(Bq)
            sgn = B >= 0.0 ? 1.0 : -1.0
            out[i] = sgn * sum(w .* (gvals .* Sval))
        end
        return out
    end

    denom = KB_MEV_PER_K * T
    for (i, B) in enumerate(Bq)
        total = 0.0
        for j in eachindex(gvals)
            arg = gvals[j] * MU_B_MEV_PER_T * B * Sval / denom
            total += w[j] * gvals[j] * Sval * tanh(arg)
        end
        out[i] = total
    end
    return out
end

# -----------------------------
# Full model
# -----------------------------

function model_magnetization(B_fields_T::AbstractVector{<:Real}, p::FitParams, draws, quad;
                             S_disp::Real = 0.5,
                             S_imp::Real = 0.5,
                             temperature_K::Real = 0.4,
                             dmax_mode::Symbol = :high_symmetry,
                             response_mode::Symbol = :linear_saturation)
    env = draw_local_env_Bsat_and_msat(p, draws; S = S_disp, dmax_mode = dmax_mode)

    Mdisp = if response_mode == :linear_saturation
        ground_magnetization_linear_saturation(B_fields_T, env)
    elseif response_mode == :smooth_crossover
        ground_magnetization_smooth_crossover(B_fields_T, env)
    else
        error("Unknown response_mode $(response_mode). Use :linear_saturation or :smooth_crossover.")
    end

    B = Float64.(B_fields_T)
    Mvv = p.chi_vv_muB_per_T .* B
    Mimp = impurity_spinhalf_magnetization(B_fields_T, p, quad;
        temperature_K = temperature_K,
        S_imp = S_imp,
    )

    Mdisp_scaled = p.scale_disp .* Mdisp
    Mvv_scaled   = p.scale_disp .* Mvv
    Mimp_scaled  = p.scale_imp  .* Mimp
    Mtotal = Mdisp_scaled .+ Mvv_scaled .+ Mimp_scaled

    return (; B_T = B,
              M_muB_per_Yb = Mtotal,
              M_total_muB_per_Yb = Mtotal,
              M_disp_scaled = Mdisp_scaled,
              M_vv_scaled = Mvv_scaled,
              M_imp_scaled = Mimp_scaled,
              M_disp_unscaled = Mdisp,
              M_vv_unscaled = Mvv,
              M_imp_unscaled = Mimp,
              env = env)
end

# -----------------------------
# Fitting
# -----------------------------

function fit_objective(x::AbstractVector{<:Real}, Bfit, Mfit, draws, quad;
                       lower, upper,
                       S_disp::Real = 0.5,
                       S_imp::Real = 0.5,
                       temperature_K::Real = 0.4,
                       dmax_mode::Symbol = :high_symmetry,
                       response_mode::Symbol = :linear_saturation)
    xcl = clamp_vector(x, lower, upper)
    p = params_from_vector(xcl)

    pred = model_magnetization(Bfit, p, draws, quad;
        S_disp = S_disp,
        S_imp = S_imp,
        temperature_K = temperature_K,
        dmax_mode = dmax_mode,
        response_mode = response_mode,
    ).M_muB_per_Yb

    if any(x -> !isfinite(x), pred)
        return 1e30
    end

    resid = pred .- Mfit
    return mean(resid .^ 2)
end

function coordinate_search(objective, x0;
                           lower,
                           upper,
                           initial_steps,
                           max_iter::Int = 180,
                           shrink::Real = 0.5,
                           min_step_factor::Real = 1e-4,
                           verbose::Bool = true)
    x = clamp_vector(x0, lower, upper)
    steps = Float64.(initial_steps)
    steps0 = copy(steps)
    f = objective(x)
    trace = NamedTuple[]

    for iter in 1:max_iter
        improved = false
        for j in eachindex(x)
            best_xj = x[j]
            best_fj = f

            for sgn in (+1.0, -1.0)
                xt = copy(x)
                xt[j] = clamp(xt[j] + sgn * steps[j], lower[j], upper[j])
                ft = objective(xt)
                if ft < best_fj
                    best_fj = ft
                    best_xj = xt[j]
                end
            end

            if best_fj < f
                x[j] = best_xj
                f = best_fj
                improved = true
            end
        end

        push!(trace, (; iter, f, x = copy(x), steps = copy(steps)))

        if verbose && (iter == 1 || iter % 10 == 0 || !improved)
            @printf("iter %4d  RMSE = %.7g  steps max = %.4g\n", iter, sqrt(f), maximum(steps))
        end

        if !improved
            steps .*= Float64(shrink)
            if maximum(steps ./ steps0) < min_step_factor
                verbose && println("Stopping: steps below tolerance.")
                break
            end
        end
    end

    return (; x, f, rmse = sqrt(f), trace)
end

# -----------------------------
# Output helpers
# -----------------------------

function write_curve_csv(filename::AbstractString, result)
    open(filename, "w") do io
        println(io, "B_T,M_total_muB_per_Yb,M_disp_scaled,M_vv_scaled,M_imp_scaled,M_disp_unscaled,M_vv_unscaled,M_imp_unscaled")
        for i in eachindex(result.B_T)
            println(io, @sprintf("%.10g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g",
                result.B_T[i],
                result.M_total_muB_per_Yb[i],
                result.M_disp_scaled[i],
                result.M_vv_scaled[i],
                result.M_imp_scaled[i],
                result.M_disp_unscaled[i],
                result.M_vv_unscaled[i],
                result.M_imp_unscaled[i]))
        end
    end
    return filename
end

function write_data_csv(filename::AbstractString, data)
    open(filename, "w") do io
        println(io, "B_T,M_muB_per_Yb")
        for i in eachindex(data.B_T)
            println(io, @sprintf("%.10g,%.12g", data.B_T[i], data.M_muB_per_Yb[i]))
        end
    end
    return filename
end

function ensemble_summary(env, p::FitParams, quad)
    finite_b = filter(isfinite, env.Bsat)
    finite_g = filter(isfinite, env.g)
    g2vals = p.gzz2 .+ p.sigma_gzz2 .* quad.z
    g2valid = g2vals[isfinite.(g2vals) .& (g2vals .> 0.0)]
    return (; mean_g = mean(finite_g),
              std_g = std(finite_g),
              median_Bsat = median(finite_b),
              q10_Bsat = quantile(finite_b, 0.10),
              q90_Bsat = quantile(finite_b, 0.90),
              mean_g2 = mean(g2valid),
              std_g2 = std(g2valid))
end

function write_fit_summary(filename::AbstractString, p0::FitParams, pbest::FitParams, fit;
                           Bmin_fit_T, Bmax_fit_T, response_mode, dmax_mode,
                           temperature_K, n_samples_fit, n_samples_final, quad_n, summary)
    open(filename, "w") do io
        println(io, "YZGO full-field magnetization fit with impurity component")
        println(io, "=========================================================")
        println(io)
        println(io, @sprintf("fit window: %.6g <= B <= %.6g T", Bmin_fit_T, Bmax_fit_T))
        println(io, @sprintf("response_mode: %s", String(response_mode)))
        println(io, @sprintf("dmax_mode: %s", String(dmax_mode)))
        println(io, @sprintf("temperature_K for impurity spin-half component: %.8g", temperature_K))
        println(io, @sprintf("n_samples_fit for dispersive component: %d", n_samples_fit))
        println(io, @sprintf("n_samples_final for dispersive component: %d", n_samples_final))
        println(io, @sprintf("normal quadrature points for impurity component: %d", quad_n))
        println(io)
        println(io, "Model")
        println(io, "-----")
        println(io, "M_total(B) = scale_disp * (M_disp(B) + chi_vv_muB_per_T * B) + scale_imp * M_imp(B)")
        println(io, "M_imp is an isolated S = 1/2 Brillouin/tanh response with mean gzz2 and sigma_gzz2.")
        println(io)
        println(io, "Initial parameters")
        println(io, "------------------")
        write_params(io, p0)
        println(io)
        println(io, "Best-fit parameters")
        println(io, "-------------------")
        write_params(io, pbest)
        println(io)
        println(io, @sprintf("RMSE in fit window: %.10g mu_B/Yb", fit.rmse))
        println(io, @sprintf("Clean Bsat from best-fit mean dispersive parameters: %.10g T", clean_Bsat_T(pbest; dmax_mode = dmax_mode)))
        println(io)
        println(io, "Final ensemble summary")
        println(io, "----------------------")
        println(io, @sprintf("dispersive mean gzz: %.10g", summary.mean_g))
        println(io, @sprintf("dispersive std gzz: %.10g", summary.std_g))
        println(io, @sprintf("Bsat median: %.10g T", summary.median_Bsat))
        println(io, @sprintf("Bsat 10-90%%: [%.10g, %.10g] T", summary.q10_Bsat, summary.q90_Bsat))
        println(io, @sprintf("impurity quadrature mean gzz2: %.10g", summary.mean_g2))
        println(io, @sprintf("impurity quadrature std gzz2: %.10g", summary.std_g2))
    end
    return filename
end

function plot_fit(experiment_all, fit_data, fit_curve;
                  Bmin_fit_T::Real = 0.0,
                  Bmax_fit_T::Real = 7.0,
                  save_png::Union{Nothing,AbstractString} = "YZGO_full_field_impurity_fit.png")
    if !HAS_GLMAKIE
        @warn "GLMakie is not available; skipping plot."
        return nothing
    end

    fig = Figure(size = (980, 620))
    ax = Axis(fig[1, 1],
        xlabel = "Magnetic field B (T)",
        ylabel = "M (mu_B / Yb)",
        title = "YZGO full-field magnetization fit with impurity component")

    scatter!(ax, experiment_all.B_T, experiment_all.M_muB_per_Yb;
        markersize = 4, label = "digitized black curve")
    scatter!(ax, fit_data.B_T, fit_data.M_muB_per_Yb;
        markersize = 6, label = @sprintf("fit data: %.1f-%.1f T", Bmin_fit_T, Bmax_fit_T))

    lines!(ax, fit_curve.B_T, fit_curve.M_total_muB_per_Yb;
        linewidth = 3, label = "total model")
    lines!(ax, fit_curve.B_T, fit_curve.M_disp_scaled;
        linewidth = 2, linestyle = :dash, label = "dispersive ground-doublet")
    lines!(ax, fit_curve.B_T, fit_curve.M_vv_scaled;
        linewidth = 2, linestyle = :dot, label = "Van Vleck term")
    lines!(ax, fit_curve.B_T, fit_curve.M_imp_scaled;
        linewidth = 2, linestyle = :dashdot, label = "non-dispersive S=1/2")

    axislegend(ax; position = :rb)
    display(fig)

    if !isnothing(save_png)
        try
            save(save_png, fig)
            println("Saved plot: $(save_png)")
        catch err
            @warn "Could not save plot" exception = (err, catch_backtrace())
        end
    end
    return fig
end

# -----------------------------
# Main driver
# -----------------------------

function fit_full_field_magnetization_impurity(; 
        experiment_csv::AbstractString = default_experiment_path(),
        Bmin_fit_T::Real = 0.0,
        Bmax_T::Real = 7.0,
        Bmax_fit_T::Real = Bmax_T,
        Bmin_plot_T::Real = 0.0,
        dB_plot_T::Real = 0.01,
        initial::FitParams = FitParams(
            scale_disp = 0.45,
            scale_imp = 0.03,
            chi_vv_muB_per_T = 0.07,
            sigma_gzz = 1.0 / 3.0,
            sigma_J1 = 1.0 / 3.0,
            gzz = 3.7,
            J1_meV = 0.177824,
            J2_meV = 0.00973612,
            gzz2 = 2.88,
            sigma_gzz2 = 0.5,
        ),
        # Bounds can be loosened if you intentionally want the optimizer to roam farther.
        lower = [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.001, 0.0, 0.2, 0.0],
        upper = [2.0, 2.0, 0.30, 2.0, 2.0, 7.0, 0.60, 0.20, 7.0, 2.5],
        initial_steps = [0.05, 0.02, 0.01, 0.05, 0.05, 0.10, 0.005, 0.001, 0.10, 0.05],
        n_samples_fit::Int = 15_000,
        n_samples_final::Int = 150_000,
        seed::Int = 20260520,
        quad_n::Int = 101,
        S_disp::Real = 0.5,
        S_imp::Real = 0.5,
        temperature_K::Real = 0.4,
        dmax_mode::Symbol = :high_symmetry,
        response_mode::Symbol = :linear_saturation,
        max_iter::Int = 180,
        verbose::Bool = true,
        write_outputs::Bool = true,
        make_plot::Bool = true)

    experiment_all = read_experimental_magnetization_csv(experiment_csv)
    fit_data = filter_fit_window(experiment_all; Bmin_fit_T = Bmin_fit_T, Bmax_fit_T = Bmax_fit_T)

    println("Loaded experimental data: $(experiment_csv)")
    println(@sprintf("Using %d points in %.6g <= B <= %.6g T", length(fit_data.B_T), Bmin_fit_T, Bmax_fit_T))
    println(@sprintf("Bmax_T for plotted curve = %.6g T", Bmax_T))
    println(@sprintf("Impurity spin-half temperature_K = %.8g", temperature_K))
    println("Initial parameters:")
    print_params(initial; prefix = "  ")
    println()

    draws_fit = make_standard_normal_draws(n_samples_fit; seed = seed)
    quad = make_normal_quadrature(quad_n)
    x0 = vector_from_params(initial)

    obj = x -> fit_objective(x, fit_data.B_T, fit_data.M_muB_per_Yb, draws_fit, quad;
        lower = lower,
        upper = upper,
        S_disp = S_disp,
        S_imp = S_imp,
        temperature_K = temperature_K,
        dmax_mode = dmax_mode,
        response_mode = response_mode,
    )

    f0 = obj(x0)
    println(@sprintf("Initial RMSE = %.8g mu_B/Yb", sqrt(f0)))
    println()

    fit = coordinate_search(obj, x0;
        lower = Float64.(lower),
        upper = Float64.(upper),
        initial_steps = Float64.(initial_steps),
        max_iter = max_iter,
        verbose = verbose,
    )

    pbest = params_from_vector(fit.x)
    println()
    println("Best-fit parameters:")
    print_params(pbest; prefix = "  ")
    println(@sprintf("Best-fit RMSE = %.8g mu_B/Yb", fit.rmse))
    println(@sprintf("Clean Bsat(best-fit mean dispersive parameters) = %.8g T", clean_Bsat_T(pbest; S = S_disp, dmax_mode = dmax_mode)))
    println()

    # Final higher-statistics curve using independent fixed draws.
    draws_final = make_standard_normal_draws(n_samples_final; seed = seed + 1)
    Bplot = collect(Float64(Bmin_plot_T):Float64(dB_plot_T):Float64(Bmax_T))
    fit_curve = model_magnetization(Bplot, pbest, draws_final, quad;
        S_disp = S_disp,
        S_imp = S_imp,
        temperature_K = temperature_K,
        dmax_mode = dmax_mode,
        response_mode = response_mode,
    )
    summary = ensemble_summary(fit_curve.env, pbest, quad)

    println("Final ensemble summary:")
    println(@sprintf("  dispersive mean gzz = %.8g, std gzz = %.8g", summary.mean_g, summary.std_g))
    println(@sprintf("  Bsat median = %.8g T, 10-90%% = [%.8g, %.8g] T",
        summary.median_Bsat, summary.q10_Bsat, summary.q90_Bsat))
    println(@sprintf("  impurity gzz2 = %.8g, sigma_gzz2 = %.8g", pbest.gzz2, pbest.sigma_gzz2))

    files = (;)
    if write_outputs
        f_curve = write_curve_csv("YZGO_full_field_impurity_fit_curve.csv", fit_curve)
        f_data = write_data_csv("YZGO_full_field_impurity_fit_data_used.csv", fit_data)
        f_summary = write_fit_summary("YZGO_full_field_impurity_fit_params.txt", initial, pbest, fit;
            Bmin_fit_T = Bmin_fit_T,
            Bmax_fit_T = Bmax_fit_T,
            response_mode = response_mode,
            dmax_mode = dmax_mode,
            temperature_K = temperature_K,
            n_samples_fit = n_samples_fit,
            n_samples_final = n_samples_final,
            quad_n = quad_n,
            summary = summary)
        println()
        println("Wrote:")
        println("  " * f_curve)
        println("  " * f_data)
        println("  " * f_summary)
        files = (; curve_csv = f_curve, fit_data_csv = f_data, summary_txt = f_summary)
    end

    fig = make_plot ? plot_fit(experiment_all, fit_data, fit_curve;
        Bmin_fit_T = Bmin_fit_T,
        Bmax_fit_T = Bmax_fit_T,
        save_png = write_outputs ? "YZGO_full_field_impurity_fit.png" : nothing) : nothing

    return (; initial,
              best = pbest,
              fit,
              experiment_all,
              fit_data,
              fit_curve,
              response_mode,
              dmax_mode,
              temperature_K = Float64(temperature_K),
              ensemble_summary = summary,
              files,
              figure = fig)
end


end # module YZGOMagnetizationFit

const NF = YZGONeutronFit
const MF = YZGOMagnetizationFit

using DelimitedFiles

using Random
using Statistics
using Printf
using Dates
using Optim
using GLMakie

const COFIT_MU_B_MEV_PER_T = 5.7883818060e-2
const COFIT_KB_MEV_PER_K   = 8.617333262e-2

# -----------------------------------------------------------------------------
# Data-loading helper missing from some monolithic versions of the neutron script
# -----------------------------------------------------------------------------

function load_neutron_fit_data_1d(; base_dir::AbstractString=NF.BASEDIR,
                                  data_mode::Symbol=:tail_bgsub,
                                  background_kwargs...)
    scans_raw = NF.load_scans(base_dir)

    if data_mode == :raw
        return scans_raw, Dict{String,Any}(), scans_raw
    elseif data_mode in (:tail_bgsub, :spline_bgsub, :bgsub)
        data_scans, background_models = NF.make_spline_background_subtracted_scans(scans_raw;
            background_kwargs...)
        return data_scans, background_models, scans_raw
    elseif data_mode == :constant_bgsub
        data_scans = NF.make_bg_subtracted_scans(scans_raw; background_kwargs...)
        return data_scans, Dict{String,Any}(), scans_raw
    else
        error("Unknown data_mode=$(data_mode). Use :tail_bgsub, :spline_bgsub, :bgsub, :constant_bgsub, or :raw.")
    end
end

# -----------------------------------------------------------------------------
# Shared parameter vector for co-fit
# -----------------------------------------------------------------------------


function load_neutron_fit_data_1d_filtered(; base_dir::AbstractString=NF.BASEDIR,
                                           data_mode::Symbol=:tail_bgsub,
                                           Ei_meV::Union{Nothing,Real}=nothing,
                                           temperature_K::Union{Nothing,Real}=nothing,
                                           atol_Ei::Real=1e-3,
                                           atol_T::Real=1e-3,
                                           background_kwargs...)
    files = filter(readdir(base_dir; join=true)) do f
        endswith(f, ".dat") && occursin("_Escan_", basename(f)) && occursin("_SYM.dat", basename(f))
    end

    scans_raw = Dict{String,Dict{Float64,NF.Scan1D}}()
    for f in sort(files)
        scan = NF.load_scan(f)
        if Ei_meV !== nothing && abs(scan.Ei_meV - Float64(Ei_meV)) > Float64(atol_Ei)
            continue
        end
        if temperature_K !== nothing && abs(scan.temperature_K - Float64(temperature_K)) > Float64(atol_T)
            continue
        end
        byfield = get!(scans_raw, scan.qtag, Dict{Float64,NF.Scan1D}())
        byfield[scan.field_T] = scan
    end

    if data_mode == :raw
        return scans_raw, Dict{String,Any}(), scans_raw
    elseif data_mode in (:tail_bgsub, :spline_bgsub, :bgsub)
        data_scans, background_models = NF.make_spline_background_subtracted_scans(scans_raw;
            background_kwargs...)
        return data_scans, background_models, scans_raw
    elseif data_mode == :constant_bgsub
        data_scans = NF.make_bg_subtracted_scans(scans_raw; background_kwargs...)
        return data_scans, Dict{String,Any}(), scans_raw
    else
        error("Unknown data_mode=$(data_mode). Use :tail_bgsub, :spline_bgsub, :bgsub, :constant_bgsub, or :raw.")
    end
end

function cofit_default_param_specs(; initial_shared_r2::Real=0.1804903296,
                                    initial_gzz::Real=3.889277345,
                                    initial_J1_meV::Real=0.2381862508,
                                    initial_J2_meV::Real=0.005661869805,
                                    initial_sigma_gzz::Real=0.3562571704,
                                    initial_sigma_J1::Real=0.1959601801,
                                    initial_sigma_J2::Real=0.3351272483,
                                    initial_gzz2::Real=3.022618403,
                                    initial_sigma_gzz2::Real=0.8370383994,
                                    initial_gperp::Union{Nothing,Real}=1.9622,
                                    initial_gperp2::Union{Nothing,Real}=5.6675,
                                    initial_chi_vv_muB_per_T::Real=0.07655025737,
                                    gzz_bounds=(1.0, 8.0),
                                    J1_bounds=(0.0, 0.60),
                                    J2_bounds=(0.0, 0.40),
                                    sigma_gzz_bounds=(0.0, 1.5),
                                    sigma_J1_bounds=(0.0, 2.0),
                                    sigma_J2_bounds=(0.0, 2.0),
                                    gzz2_bounds=(0.5, 8.0),
                                    sigma_gzz2_bounds=(0.0, 2.0),
                                    gperp_bounds=(0.0, 8.0),
                                    gperp2_bounds=(0.0, 8.0),
                                    chi_vv_bounds=(0.0, 0.30),
                                    log10_shared_r2_bounds=(-4.0, 4.0))
    pars = NF.demo_defaults()
    return NF.FitParamSpec[
        NF.FitParamSpec(:gzz,       Float64(gzz_bounds[1]),       Float64(gzz_bounds[2]),       Float64(initial_gzz)),
        NF.FitParamSpec(:J1_meV,    Float64(J1_bounds[1]),        Float64(J1_bounds[2]),        Float64(initial_J1_meV)),
        NF.FitParamSpec(:J2_meV,    Float64(J2_bounds[1]),        Float64(J2_bounds[2]),        Float64(initial_J2_meV)),
        NF.FitParamSpec(:sigma_gzz, Float64(sigma_gzz_bounds[1]), Float64(sigma_gzz_bounds[2]), Float64(initial_sigma_gzz)),
        NF.FitParamSpec(:sigma_J1,  Float64(sigma_J1_bounds[1]),  Float64(sigma_J1_bounds[2]),  Float64(initial_sigma_J1)),
        NF.FitParamSpec(:sigma_J2,  Float64(sigma_J2_bounds[1]),  Float64(sigma_J2_bounds[2]),  Float64(initial_sigma_J2)),
        NF.FitParamSpec(:gzz2,      Float64(gzz2_bounds[1]),      Float64(gzz2_bounds[2]),      Float64(initial_gzz2)),
        NF.FitParamSpec(:sigma_gzz2,Float64(sigma_gzz2_bounds[1]),Float64(sigma_gzz2_bounds[2]),Float64(initial_sigma_gzz2)),
        NF.FitParamSpec(:gperp,     Float64(gperp_bounds[1]),     Float64(gperp_bounds[2]),     Float64(initial_gperp === nothing ? pars.model.gzz : initial_gperp)),
        NF.FitParamSpec(:gperp2,    Float64(gperp2_bounds[1]),    Float64(gperp2_bounds[2]),    Float64(initial_gperp2 === nothing ? initial_gzz2 : initial_gperp2)),
        NF.FitParamSpec(:chi_vv_muB_per_T,
                        Float64(chi_vv_bounds[1]), Float64(chi_vv_bounds[2]),
                        Float64(initial_chi_vv_muB_per_T)),
        NF.FitParamSpec(:log10_second_kernel_relative_intensity,
                        Float64(log10_shared_r2_bounds[1]),
                        Float64(log10_shared_r2_bounds[2]),
                        log10(Float64(initial_shared_r2))),
    ]
end

function neutron_param_dict_from_cofit(p::Dict{Symbol,Float64})
    pn = Dict{Symbol,Float64}()
    for nm in [:gzz, :J1_meV, :J2_meV, :sigma_gzz, :sigma_J1, :sigma_J2, :gzz2, :sigma_gzz2]
        pn[nm] = p[nm]
    end
    if haskey(p, :gperp)
        pn[:gperp] = p[:gperp]
    end
    if haskey(p, :gperp2)
        pn[:gperp2] = p[:gperp2]
    end
    pn[:log10_second_kernel_relative_intensity] = p[:log10_second_kernel_relative_intensity]
    return pn
end

shared_r2_from_cofit(p::Dict{Symbol,Float64}) = 10.0 ^ p[:log10_second_kernel_relative_intensity]
neutron_r2_from_cofit(p::Dict{Symbol,Float64}) = shared_r2_from_cofit(p)
magnetization_r2_from_cofit(p::Dict{Symbol,Float64}) = shared_r2_from_cofit(p)

# -----------------------------------------------------------------------------
# Magnetization model with sigma_J2 no longer forced equal to sigma_J1
# -----------------------------------------------------------------------------

function make_standard_normal_draws_cofit(n_samples::Int; seed::Int=20260520)
    rng = MersenneTwister(seed)
    return (; zg = randn(rng, n_samples), z1 = randn(rng, n_samples), z2 = randn(rng, n_samples))
end

function make_normal_quadrature_cofit(n::Int=101; zmax::Real=5.0)
    n >= 5 || error("normal quadrature needs at least 5 points")
    z = collect(range(-Float64(zmax), Float64(zmax); length=n))
    w = exp.(-0.5 .* z .^ 2)
    w ./= sum(w)
    return (; z, w)
end

function exchange_Dmax_meV_cofit(J1_meV::Real, J2_meV::Real; mode::Symbol=:high_symmetry)
    J1 = Float64(J1_meV)
    J2 = Float64(J2_meV)
    if mode == :high_symmetry
        return max(0.0, 9.0 * J1, 8.0 * (J1 + J2))
    elseif mode == :K_only
        return max(0.0, 9.0 * J1)
    elseif mode == :M_only
        return max(0.0, 8.0 * (J1 + J2))
    else
        error("Unknown dmax mode $(mode). Use :high_symmetry, :K_only, or :M_only.")
    end
end

function draw_local_env_Bsat_and_msat_cofit(p::Dict{Symbol,Float64}, draws;
                                            S::Real=0.5,
                                            dmax_mode::Symbol=:high_symmetry)
    n = length(draws.zg)
    g = Vector{Float64}(undef, n)
    Bsat = Vector{Float64}(undef, n)
    msat = Vector{Float64}(undef, n)

    sJ1_abs = p[:sigma_J1] * abs(p[:J1_meV])
    sJ2_abs = p[:sigma_J2] * abs(p[:J2_meV])
    Sval = Float64(S)

    for i in 1:n
        gi = p[:gzz] + p[:sigma_gzz] * draws.zg[i]
        J1i = p[:J1_meV] + sJ1_abs * draws.z1[i]
        J2i = p[:J2_meV] + sJ2_abs * draws.z2[i]
        Dmaxi = exchange_Dmax_meV_cofit(J1i, J2i; mode=dmax_mode)

        g[i] = gi
        msat[i] = gi * Sval
        if isfinite(gi) && gi > 0.0 && isfinite(Dmaxi) && Dmaxi > 0.0
            Bsat[i] = Sval * Dmaxi / (gi * COFIT_MU_B_MEV_PER_T)
        elseif isfinite(gi) && gi > 0.0
            Bsat[i] = 0.0
        else
            Bsat[i] = NaN
        end
    end
    return (; g, Bsat, msat)
end

function ground_magnetization_linear_saturation_cofit(B_fields_T::AbstractVector{<:Real}, env)
    Bq = Float64.(B_fields_T)
    valid = isfinite.(env.Bsat) .& isfinite.(env.msat) .& (env.g .> 0.0)
    nvalid = count(valid)
    nvalid == 0 && return fill(NaN, length(Bq))

    b = env.Bsat[valid]
    m = env.msat[valid]

    zero_or_negative = b .<= 0.0
    base_sat = sum(m[zero_or_negative])

    bpos = b[.!zero_or_negative]
    mpos = m[.!zero_or_negative]

    if isempty(bpos)
        return fill(base_sat / nvalid, length(Bq))
    end

    idx = sortperm(bpos)
    bs = bpos[idx]
    ms = mpos[idx]
    prefix_m = cumsum(ms)
    slopes = ms ./ bs
    suffix_slope = reverse(cumsum(reverse(slopes)))

    out = similar(Bq)
    for (i, B) in enumerate(Bq)
        sgn = B >= 0.0 ? 1.0 : -1.0
        Ba = abs(B)
        k = searchsortedlast(bs, Ba)
        saturated_sum = k > 0 ? prefix_m[k] : 0.0
        unsat_slope_sum = k < length(bs) ? suffix_slope[k + 1] : 0.0
        out[i] = sgn * (base_sat + saturated_sum + Ba * unsat_slope_sum) / nvalid
    end
    return out
end

function ground_magnetization_smooth_crossover_cofit(B_fields_T::AbstractVector{<:Real}, env)
    Bq = Float64.(B_fields_T)
    valid = isfinite.(env.Bsat) .& isfinite.(env.msat) .& (env.g .> 0.0)
    nvalid = count(valid)
    nvalid == 0 && return fill(NaN, length(Bq))

    b = env.Bsat[valid]
    m = env.msat[valid]
    out = zeros(Float64, length(Bq))
    for (i, B) in enumerate(Bq)
        total = 0.0
        for j in eachindex(b)
            if b[j] <= 0.0
                total += B >= 0.0 ? m[j] : -m[j]
            else
                total += m[j] * B / sqrt(B^2 + b[j]^2)
            end
        end
        out[i] = total / nvalid
    end
    return out
end

function nondispersive_spinhalf_magnetization_cofit(B_fields_T::AbstractVector{<:Real}, p::Dict{Symbol,Float64}, quad;
                                                    temperature_K::Real=0.4,
                                                    S_imp::Real=0.5)
    Bq = Float64.(B_fields_T)
    T = Float64(temperature_K)
    Sval = Float64(S_imp)

    gvals_all = p[:gzz2] .+ p[:sigma_gzz2] .* quad.z
    valid = isfinite.(gvals_all) .& (gvals_all .> 0.0)
    count(valid) == 0 && return fill(NaN, length(Bq))

    gvals = gvals_all[valid]
    w = quad.w[valid]
    w ./= sum(w)

    out = zeros(Float64, length(Bq))
    if T <= 0.0
        for (i, B) in enumerate(Bq)
            sgn = B >= 0.0 ? 1.0 : -1.0
            out[i] = sgn * sum(w .* (gvals .* Sval))
        end
        return out
    end

    denom = COFIT_KB_MEV_PER_K * T
    for (i, B) in enumerate(Bq)
        total = 0.0
        for j in eachindex(gvals)
            arg = gvals[j] * COFIT_MU_B_MEV_PER_T * B * Sval / denom
            total += w[j] * gvals[j] * Sval * tanh(arg)
        end
        out[i] = total
    end
    return out
end

function magnetization_components_cofit(B_fields_T::AbstractVector{<:Real}, p::Dict{Symbol,Float64}, draws, quad;
                                        S_disp::Real=0.5,
                                        S_imp::Real=0.5,
                                        temperature_K::Real=0.4,
                                        dmax_mode::Symbol=:high_symmetry,
                                        response_mode::Symbol=:linear_saturation)
    env = draw_local_env_Bsat_and_msat_cofit(p, draws; S=S_disp, dmax_mode=dmax_mode)

    Mdisp = if response_mode == :linear_saturation
        ground_magnetization_linear_saturation_cofit(B_fields_T, env)
    elseif response_mode == :smooth_crossover
        ground_magnetization_smooth_crossover_cofit(B_fields_T, env)
    else
        error("Unknown response_mode $(response_mode). Use :linear_saturation or :smooth_crossover.")
    end

    B = Float64.(B_fields_T)
    Mvv = p[:chi_vv_muB_per_T] .* B
    Mnon = nondispersive_spinhalf_magnetization_cofit(B_fields_T, p, quad;
        temperature_K=temperature_K,
        S_imp=S_imp,
    )
    r2 = magnetization_r2_from_cofit(p)
    Mcombo = Mdisp .+ Mvv .+ r2 .* Mnon

    return (; B_T=B,
              M_combo_unscaled=Mcombo,
              M_disp_unscaled=Mdisp,
              M_vv_unscaled=Mvv,
              M_nondispersive_unscaled=Mnon,
              r2=r2,
              env=env)
end

function best_positive_scale(y::AbstractVector{<:Real}, m::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    denom = sum(w .* Float64.(m) .* Float64.(m))
    denom > 0 || return NaN
    s = sum(w .* Float64.(y) .* Float64.(m)) / denom
    return max(s, 0.0)
end

function magnetization_sigma_vector(Mfit::AbstractVector{<:Real};
                                    magnetization_sigma_muB::Union{Nothing,Real}=nothing,
                                    magnetization_error_fraction_of_range::Real=0.02,
                                    magnetization_error_floor_muB::Real=1e-4)
    Mf = Float64.(Mfit)
    if magnetization_sigma_muB !== nothing
        sigma = Float64(magnetization_sigma_muB)
    else
        rng = maximum(Mf) - minimum(Mf)
        sigma = Float64(magnetization_error_fraction_of_range) * rng
    end
    sigma = max(sigma, Float64(magnetization_error_floor_muB), eps(Float64))
    return fill(sigma, length(Mf))
end

function evaluate_magnetization_fit(p::Dict{Symbol,Float64}, mag_data, draws, quad;
                                    S_disp::Real=0.5,
                                    S_imp::Real=0.5,
                                    temperature_K::Real=0.4,
                                    dmax_mode::Symbol=:high_symmetry,
                                    response_mode::Symbol=:linear_saturation,
                                    magnetization_sigma_muB::Union{Nothing,Real}=nothing,
                                    magnetization_error_fraction_of_range::Real=0.02,
                                    magnetization_error_floor_muB::Real=1e-4,
                                    nfree::Int=1,
                                    kwargs...)
    comp = magnetization_components_cofit(mag_data.B_T, p, draws, quad;
        S_disp=S_disp,
        S_imp=S_imp,
        temperature_K=temperature_K,
        dmax_mode=dmax_mode,
        response_mode=response_mode,
    )

    if any(x -> !isfinite(x), comp.M_combo_unscaled)
        return (; redchi2=Inf, scale=NaN, components=comp, sigma=Float64[])
    end

    sigma = magnetization_sigma_vector(mag_data.M_muB_per_Yb;
        magnetization_sigma_muB=magnetization_sigma_muB,
        magnetization_error_fraction_of_range=magnetization_error_fraction_of_range,
        magnetization_error_floor_muB=magnetization_error_floor_muB,
    )
    w = 1.0 ./ sigma.^2
    scale = best_positive_scale(mag_data.M_muB_per_Yb, comp.M_combo_unscaled, w)
    if !isfinite(scale) || scale < 0
        return (; redchi2=Inf, scale=scale, components=comp, sigma=sigma)
    end

    resid = Float64.(mag_data.M_muB_per_Yb) .- scale .* comp.M_combo_unscaled
    chi2 = sum(w .* resid.^2)
    dof = max(length(resid) - nfree, 1)
    return (; redchi2=chi2 / dof, scale=scale, components=comp, sigma=sigma)
end

function scaled_magnetization_curve(p::Dict{Symbol,Float64}, Bgrid::AbstractVector{<:Real}, scale::Real, draws, quad;
                                    S_disp::Real=0.5,
                                    S_imp::Real=0.5,
                                    temperature_K::Real=0.4,
                                    dmax_mode::Symbol=:high_symmetry,
                                    response_mode::Symbol=:linear_saturation)
    comp = magnetization_components_cofit(Bgrid, p, draws, quad;
        S_disp=S_disp,
        S_imp=S_imp,
        temperature_K=temperature_K,
        dmax_mode=dmax_mode,
        response_mode=response_mode,
    )
    s = Float64(scale)
    return (; B_T=comp.B_T,
              M_total_muB_per_Yb=s .* comp.M_combo_unscaled,
              M_disp_scaled=s .* comp.M_disp_unscaled,
              M_vv_scaled=s .* comp.M_vv_unscaled,
              M_nondispersive_scaled=s .* comp.r2 .* comp.M_nondispersive_unscaled,
              M_disp_unscaled=comp.M_disp_unscaled,
              M_vv_unscaled=comp.M_vv_unscaled,
              M_nondispersive_unscaled=comp.M_nondispersive_unscaled,
              scale=s,
              r2=comp.r2,
              env=comp.env)
end

# -----------------------------------------------------------------------------
# Neutron objective helper
# -----------------------------------------------------------------------------

function evaluate_neutron_fit(p::Dict{Symbol,Float64}, data_scans;
                              fields_T::AbstractVector{<:Real}=[9.0],
                              qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                              fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
                              n_samples_per_cut::Int=100_000,
                              seed::Int=2026,
                              use_errors::Bool=true,
                              error_floor::Real=0.0,
                              lattice::NF.LatticeParams=NF.demo_defaults().lattice,
                              resolution::NF.ResolutionParams=NF.demo_defaults().resolution,
                              cuts::Vector{NF.CutSpec1D}=NF.default_cuts_1d(),
                              S::Real=0.5,
                              J_units::Symbol=:fractional,
                              correlate_J1_J2::Bool=false,
                              use_form_factor::Bool=true,
                              include_j2_formfactor::Bool=true,
                              include_kfki::Bool=true,
                              polarization::Symbol=:transverse_c,
                              nfree::Int=1,
                              kwargs...)
    pn = neutron_param_dict_from_cofit(p)
    r2n = NF.second_kernel_relative_intensity(pn)

    model = NF.simulate_two_kernel_model_fields(pn;
        fields_T=fields_T,
        cuts=cuts,
        lattice=lattice,
        resolution=resolution,
        n_samples_per_cut=n_samples_per_cut,
        seed=seed,
        S=S,
        J_units=J_units,
        correlate_J1_J2=correlate_J1_J2,
        use_form_factor=use_form_factor,
        include_j2_formfactor=include_j2_formfactor,
        include_kfki=include_kfki,
        polarization=polarization,
    )

    bundle = NF.build_fit_vector_bundle(data_scans, model;
        fields_T=fields_T,
        qtags=qtags,
        fit_windows_by_q=fit_windows_by_q,
        use_errors=use_errors,
        error_floor=error_floor,
    )

    if length(bundle.data_intensity) < 6
        return (; redchi2=Inf, scale=NaN, r2=r2n, model=model, bundle=bundle)
    end

    scale = NF.best_global_scale(bundle, r2n)
    if !isfinite(scale) || scale < 0
        return (; redchi2=Inf, scale=scale, r2=r2n, model=model, bundle=bundle)
    end

    redchi2 = NF.chisq_for_bundle(bundle, scale, r2n; nfree=nfree, reduced=true)
    return (; redchi2=redchi2, scale=scale, r2=r2n, model=model, bundle=bundle)
end

function evaluate_cofit_objective(u::AbstractVector{<:Real}, specs::Vector{NF.FitParamSpec},
                                  data_scans, mag_data, mag_draws, mag_quad;
                                  neutron_weight::Real=1.0,
                                  magnetization_weight::Real=1.0,
                                  penalty::Real=1e30,
                                  kwargs...)
    try
        p = NF.unpack_unconstrained(u, specs)

        nres = evaluate_neutron_fit(p, data_scans; kwargs...)
        mres = evaluate_magnetization_fit(p, mag_data, mag_draws, mag_quad; kwargs...)

        if !isfinite(nres.redchi2) || !isfinite(mres.redchi2)
            return (; total=penalty, p=p, neutron=nres, magnetization=mres)
        end

        total = Float64(neutron_weight) * nres.redchi2 + Float64(magnetization_weight) * mres.redchi2
        total = isfinite(total) ? total : penalty
        return (; total=total, p=p, neutron=nres, magnetization=mres)
    catch err
        @warn "Co-fit objective evaluation failed" exception=(err, catch_backtrace())
        p = Dict{Symbol,Float64}()
        dummy = (; redchi2=Inf, scale=NaN)
        return (; total=penalty, p=p, neutron=dummy, magnetization=dummy)
    end
end

# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------

function write_cofit_summary(path::AbstractString;
                             specs::Vector{NF.FitParamSpec},
                             params::Dict{Symbol,Float64},
                             neutron_weight::Real,
                             magnetization_weight::Real,
                             neutron_result,
                             magnetization_result,
                             fields_T,
                             qtags,
                             fit_windows_by_q,
                             Bmin_fit_T,
                             Bmax_fit_T,
                             magnetization_sigma_muB,
                             magnetization_error_fraction_of_range,
                             n_samples_per_cut,
                             final_n_samples_per_cut,
                             n_samples_magnetization,
                             seed,
                             mag_seed,
                             optimizer_result=nothing)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "YZGO neutron + magnetization co-fit summary")
        println(io, "Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io)
        println(io, "Model:")
        println(io, "  Shared magnetization/energy: gzz, J1_meV, J2_meV, sigma_gzz, sigma_J1, sigma_J2, gzz2, sigma_gzz2, r2")
        println(io, "  Neutron-only intensity prefactors: I_disp ∝ gperp^2, I_nondispersive ∝ gperp2^2")
        println(io, "  Neutron:       I = scale_neutron * (gperp^2 * I_disp + r2 * gperp2^2 * I_nondispersive)")
        println(io, "  Magnetization: M = scale_mag * (M_disp + chi_vv*B + r2 * M_nondispersive)")
        println(io, "  The same r2 is used for neutron and magnetization. This is the deliberately constrained shared-fraction variant.")
        println(io)
        println(io, @sprintf("Objective = %.8g * neutron_redchi2 + %.8g * magnetization_redchi2",
            Float64(neutron_weight), Float64(magnetization_weight)))
        println(io)
        println(io, "Neutron fields fitted (T): ", collect(fields_T))
        println(io, "Note: in this version the default fitted neutron fields are [9.0, 14.0], so the 14 T Ei=4.65 meV data are included in the objective, not only plotted as a prediction.")
        println(io, "Neutron Q cuts fitted: ", collect(qtags))
        println(io, "Neutron fit windows by qtag:")
        for qtag in qtags
            println(io, "  ", qtag, " => ", fit_windows_by_q === nothing ? nothing : get(fit_windows_by_q, qtag, nothing))
        end
        println(io)
        println(io, @sprintf("Magnetization fit window: %.8g <= B <= %.8g T", Float64(Bmin_fit_T), Float64(Bmax_fit_T)))
        if magnetization_sigma_muB === nothing
            println(io, @sprintf("Magnetization sigma: %.8g * data range", Float64(magnetization_error_fraction_of_range)))
        else
            println(io, @sprintf("Magnetization sigma: %.8g mu_B/Yb", Float64(magnetization_sigma_muB)))
        end
        println(io)
        println(io, "Monte Carlo samples per cut during optimization: ", n_samples_per_cut)
        println(io, "Monte Carlo samples per cut for final neutron model: ", final_n_samples_per_cut)
        println(io, "Magnetization disorder samples: ", n_samples_magnetization)
        println(io, "Neutron common-random-number seed: ", seed)
        println(io, "Magnetization common-random-number seed: ", mag_seed)
        println(io)
        println(io, "Parameter bounds and initials:")
        for s in specs
            println(io, @sprintf("  %-56s initial=% .10g bounds=[% .10g, % .10g]",
                                 String(s.name), s.initial, s.lo, s.hi))
        end
        println(io)
        println(io, "Best-fit parameters:")
        for nm in [:gzz, :J1_meV, :J2_meV, :sigma_gzz, :sigma_J1, :sigma_J2,
                   :gzz2, :sigma_gzz2, :gperp, :gperp2, :chi_vv_muB_per_T,
                   :log10_second_kernel_relative_intensity]
            println(io, @sprintf("  %-56s % .12g", String(nm), params[nm]))
        end
        println(io, @sprintf("  %-56s % .12g", "r2_shared", shared_r2_from_cofit(params)))
        println(io, @sprintf("  %-56s % .12g", "neutron_global_scale", neutron_result.scale))
        println(io, @sprintf("  %-56s % .12g", "magnetization_global_scale", magnetization_result.scale))
        println(io, @sprintf("  %-56s % .12g", "neutron_nondispersive_effective_scale", neutron_result.scale * neutron_result.r2))
        println(io, @sprintf("  %-56s % .12g", "magnetization_nondispersive_effective_scale", magnetization_result.scale * magnetization_result.components.r2))
        println(io)
        println(io, @sprintf("Neutron reduced chi^2: %.12g", neutron_result.redchi2))
        println(io, @sprintf("Magnetization reduced chi^2: %.12g", magnetization_result.redchi2))
        println(io, @sprintf("Weighted objective: %.12g", Float64(neutron_weight) * neutron_result.redchi2 + Float64(magnetization_weight) * magnetization_result.redchi2))
        if optimizer_result !== nothing
            println(io)
            println(io, "Optimizer summary:")
            println(io, optimizer_result)
        end
    end
    return path
end

function write_magnetization_model_csv(path::AbstractString, curve)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "B_T,M_total_muB_per_Yb,M_disp_scaled,M_vv_scaled,M_nondispersive_scaled,M_disp_unscaled,M_vv_unscaled,M_nondispersive_unscaled,magnetization_scale,r2_shared")
        for i in eachindex(curve.B_T)
            println(io, join((
                curve.B_T[i],
                curve.M_total_muB_per_Yb[i],
                curve.M_disp_scaled[i],
                curve.M_vv_scaled[i],
                curve.M_nondispersive_scaled[i],
                curve.M_disp_unscaled[i],
                curve.M_vv_unscaled[i],
                curve.M_nondispersive_unscaled[i],
                curve.scale,
                curve.r2,
            ), ","))
        end
    end
    return path
end

function plot_yzgo_cofit_combo(data_scans, model_results_by_field, experiment_all, fit_data_mag, mag_curve;
                               neutron_scale::Real,
                               neutron_r2::Real,
                               fields_T::AbstractVector{<:Real}=[9.0],
                               qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                               qtag_to_modelkey=NF.OVERLAY_QTAG_TO_MODELKEY,
                               fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
                               neutron_ylims::Union{Nothing,Tuple{Float64,Float64}}=NF.YZGO_FIT_PLOT_YLIMS,
                               extra_data_scans=nothing,
                               extra_model_results_by_field=nothing,
                               extra_fields_T::AbstractVector{<:Real}=Float64[],
                               extra_neutron_scale::Union{Nothing,Real}=nothing,
                               extra_neutron_r2::Union{Nothing,Real}=nothing,
                               extra_row_label::AbstractString="Ei = 4.65 meV prediction",
                               figure_title::AbstractString="YZGO co-fit: neutron scans + magnetization",
                               outpath::Union{Nothing,AbstractString}=nothing,
                               save_png::Bool=true,
                               display_fig::Bool=true)
    primary_fields = Float64.(collect(fields_T))
    extra_fields = Float64.(collect(extra_fields_T))
    nrows_primary = max(length(primary_fields), 1)
    nrows_extra = (extra_data_scans === nothing || extra_model_results_by_field === nothing) ? 0 : length(extra_fields)
    nrows = nrows_primary + nrows_extra
    ncols = length(qtags) + 1
    fig = Figure(size=(430*ncols, 360*nrows + 110))
    Label(fig[0, :], figure_title; fontsize=20)

    function plot_neutron_row!(irow::Int, Bf::Float64, row_data_scans, row_model_results_by_field;
                               row_scale::Real, row_r2::Real, prefix::AbstractString="")
        haskey(row_model_results_by_field, Bf) || return
        model_components = row_model_results_by_field[Bf]
        model_disp_results = model_components.dispersive
        model_flat_results = model_components.nondispersive

        for (icol, qtag) in enumerate(qtags)
            if !(haskey(row_data_scans, qtag) && haskey(row_data_scans[qtag], Bf))
                continue
            end
            scan = row_data_scans[qtag][Bf]
            title = get(NF.OVERLAY_QTAG_TITLES, qtag, qtag)
            row_title = isempty(prefix) ? @sprintf("%.0f T, %s", Bf, title) : @sprintf("%s, %.0f T, %s", prefix, Bf, title)
            ax = Axis(fig[irow, icol];
                title=row_title,
                xlabel="Delta E (meV)",
                ylabel=icol == 1 ? "Intensity - background" : "",
            )

            scatter!(ax, scan.energy, scan.intensity; label="data")
            errorbars!(ax, scan.energy, scan.intensity, scan.error; whiskerwidth=4)

            mres_disp = NF._model_result_for_qtag(model_disp_results, qtag; qtag_to_modelkey=qtag_to_modelkey)
            mres_flat = NF._model_result_for_qtag(model_flat_results, qtag; qtag_to_modelkey=qtag_to_modelkey)
            E = mres_disp.E_centers_meV
            y_disp = Float64(row_scale) .* mres_disp.intensity
            y_flat = Float64(row_scale) .* Float64(row_r2) .* mres_flat.intensity
            y_total = y_disp .+ y_flat

            lines!(ax, E, y_disp; label="dispersive", linewidth=2, linestyle=:dash)
            lines!(ax, E, y_flat; label="non-dispersive", linewidth=2, linestyle=:dot)
            lines!(ax, E, y_total; label="total", linewidth=3)

            windows = fit_windows_by_q === nothing ? Tuple{Float64,Float64}[] : get(fit_windows_by_q, qtag, Tuple{Float64,Float64}[])
            for (lo, hi) in windows
                vlines!(ax, [lo, hi]; linestyle=:dash, linewidth=1)
            end

            if hasproperty(scan, :meta)
                xlims!(ax, scan.meta.de_range...)
            end
            if neutron_ylims !== nothing
                ylims!(ax, neutron_ylims...)
            end
            axislegend(ax; position=:rt, framevisible=false)
        end
    end

    for (irow, Bf) in enumerate(primary_fields)
        plot_neutron_row!(irow, Bf, data_scans, model_results_by_field;
            row_scale=neutron_scale,
            row_r2=neutron_r2,
            prefix="")
    end

    if nrows_extra > 0
        row_scale = extra_neutron_scale === nothing ? neutron_scale : extra_neutron_scale
        row_r2 = extra_neutron_r2 === nothing ? neutron_r2 : extra_neutron_r2
        for (jrow, Bf) in enumerate(extra_fields)
            plot_neutron_row!(nrows_primary + jrow, Bf, extra_data_scans, extra_model_results_by_field;
                row_scale=row_scale,
                row_r2=row_r2,
                prefix=extra_row_label)
        end
    end

    axm = Axis(fig[1:nrows, ncols];
        title="Magnetization",
        xlabel="Magnetic field B (T)",
        ylabel="M (mu_B / Yb)",
    )
    scatter!(axm, experiment_all.B_T, experiment_all.M_muB_per_Yb;
        markersize=4, label="digitized data")
    scatter!(axm, fit_data_mag.B_T, fit_data_mag.M_muB_per_Yb;
        markersize=6, label="fit data")
    lines!(axm, mag_curve.B_T, mag_curve.M_total_muB_per_Yb;
        linewidth=3, label="total model")
    lines!(axm, mag_curve.B_T, mag_curve.M_disp_scaled;
        linewidth=2, linestyle=:dash, label="dispersive")
    lines!(axm, mag_curve.B_T, mag_curve.M_vv_scaled;
        linewidth=2, linestyle=:dot, label="Van Vleck")
    lines!(axm, mag_curve.B_T, mag_curve.M_nondispersive_scaled;
        linewidth=2, linestyle=:dashdot, label="non-dispersive S=1/2")
    axislegend(axm; position=:rb, framevisible=false)

    if save_png && outpath !== nothing
        mkpath(dirname(outpath))
        save(outpath, fig)
        println("Saved co-fit combo figure to: ", outpath)
    end
    if display_fig
        display(fig)
    end
    return fig
end

# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------

function run_yzgo_neutron_magnetization_cofit_shared_fraction(; 
        base_dir::AbstractString=NF.BASEDIR,
        magnetization_csv::AbstractString=MF.default_experiment_path(),
        outdir::AbstractString=joinpath(base_dir, "cofit_neutron_magnetization_shared_fraction_fit9T14T_4p65"),
        fields_T::AbstractVector{<:Real}=[9.0, 14.0],
        neutron_fit_Ei_meV::Union{Nothing,Real}=4.65,
        neutron_fit_temperature_K::Union{Nothing,Real}=0.07,
        qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
        fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
        data_mode::Symbol=:tail_bgsub,
        neutron_weight::Real=1.0,
        magnetization_weight::Real=10.0,
        specs::Union{Nothing,Vector{NF.FitParamSpec}}=nothing,
        maxiters::Int=1000,
        n_samples_per_cut::Int=100_000,
        final_n_samples_per_cut::Int=500_000,
        seed::Int=2026,
        n_samples_magnetization::Int=15_000,
        final_n_samples_magnetization::Int=150_000,
        mag_seed::Int=20260520,
        quad_n::Int=101,
        Bmin_fit_T::Real=0.0,
        Bmax_fit_T::Real=7.0,
        Bmin_plot_T::Real=0.0,
        Bmax_plot_T::Real=7.0,
        dB_plot_T::Real=0.01,
        magnetization_sigma_muB::Union{Nothing,Real}=nothing,
        magnetization_error_fraction_of_range::Real=0.02,
        magnetization_error_floor_muB::Real=1e-4,
        use_errors::Bool=true,
        error_floor::Real=0.0,
        show_trace::Bool=true,
        trace_every::Int=5,
        make_plots::Bool=true,
        display_figures::Bool=true,
        run_optimization::Bool=true,
        add_14T_4p65_comparison_row::Bool=false,
        comparison_Ei_meV::Real=4.65,
        comparison_temperature_K::Real=0.07,
        comparison_fields_T::AbstractVector{<:Real}=[14.0],
        neutron_plot_ylims::Union{Nothing,Tuple{Float64,Float64}}=NF.YZGO_FIT_PLOT_YLIMS,
        lattice::NF.LatticeParams=NF.demo_defaults().lattice,
        resolution::NF.ResolutionParams=NF.demo_defaults().resolution,
        cuts::Vector{NF.CutSpec1D}=NF.default_cuts_1d(),
        S::Real=0.5,
        S_imp::Real=0.5,
        temperature_K::Real=0.4,
        J_units::Symbol=:fractional,
        correlate_J1_J2::Bool=false,
        use_form_factor::Bool=true,
        include_j2_formfactor::Bool=true,
        include_kfki::Bool=true,
        polarization::Symbol=:transverse_c,
        dmax_mode::Symbol=:high_symmetry,
        response_mode::Symbol=:linear_saturation,
        background_kwargs...)

    mkpath(outdir)

    specs2 = specs === nothing ? cofit_default_param_specs() : specs
    NF.print_fit_param_specs(specs2)

    println()
    println("Loading neutron experimental scans from: ", base_dir)
    if neutron_fit_Ei_meV === nothing && neutron_fit_temperature_K === nothing
        data_scans, background_models, scans_raw = load_neutron_fit_data_1d(;
            base_dir=base_dir,
            data_mode=data_mode,
            background_kwargs...
        )
    else
        println(@sprintf("Filtering fitted neutron scans to Ei=%.5g meV, T=%.5g K",
            Float64(neutron_fit_Ei_meV), Float64(neutron_fit_temperature_K)))
        data_scans, background_models, scans_raw = load_neutron_fit_data_1d_filtered(;
            base_dir=base_dir,
            data_mode=data_mode,
            Ei_meV=neutron_fit_Ei_meV,
            temperature_K=neutron_fit_temperature_K,
            background_kwargs...
        )
    end

    println("Loading magnetization data from: ", magnetization_csv)
    mag_all = MF.read_experimental_magnetization_csv(magnetization_csv)
    mag_fit = MF.filter_fit_window(mag_all; Bmin_fit_T=Bmin_fit_T, Bmax_fit_T=Bmax_fit_T)

    println("Using neutron data mode: ", data_mode)
    println("Fitting neutron fields: ", collect(fields_T))
    println("Fitting neutron q cuts: ", collect(qtags))
    if neutron_fit_Ei_meV !== nothing || neutron_fit_temperature_K !== nothing
        println("Fitting neutron scan filter: Ei=", neutron_fit_Ei_meV, " meV, T=", neutron_fit_temperature_K, " K")
    end
    println("Fitting magnetization window: ", Bmin_fit_T, " to ", Bmax_fit_T, " T")
    println(@sprintf("Objective weights: neutron=%.8g, magnetization=%.8g", Float64(neutron_weight), Float64(magnetization_weight)))
    if magnetization_sigma_muB === nothing
        println(@sprintf("Magnetization uncertainty model: sigma = %.6g * fitted M range, floor %.6g mu_B/Yb",
            Float64(magnetization_error_fraction_of_range), Float64(magnetization_error_floor_muB)))
    else
        println(@sprintf("Magnetization uncertainty model: sigma = %.6g mu_B/Yb", Float64(magnetization_sigma_muB)))
    end
    println()

    u0 = NF.unconstrained_initial(specs2)
    mag_draws = make_standard_normal_draws_cofit(n_samples_magnetization; seed=mag_seed)
    mag_quad = make_normal_quadrature_cofit(quad_n)

    eval_counter = Ref(0)
    last_metrics = Ref{Any}(nothing)

    obj = function(u)
        eval_counter[] += 1
        metrics = evaluate_cofit_objective(u, specs2, data_scans, mag_fit, mag_draws, mag_quad;
            neutron_weight=neutron_weight,
            magnetization_weight=magnetization_weight,
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            n_samples_per_cut=n_samples_per_cut,
            seed=seed,
            use_errors=use_errors,
            error_floor=error_floor,
            lattice=lattice,
            resolution=resolution,
            cuts=cuts,
            S=S,
            S_disp=S,
            S_imp=S_imp,
            temperature_K=temperature_K,
            J_units=J_units,
            correlate_J1_J2=correlate_J1_J2,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
            dmax_mode=dmax_mode,
            response_mode=response_mode,
            magnetization_sigma_muB=magnetization_sigma_muB,
            magnetization_error_fraction_of_range=magnetization_error_fraction_of_range,
            magnetization_error_floor_muB=magnetization_error_floor_muB,
        )
        last_metrics[] = metrics

        if trace_every > 0 && (eval_counter[] == 1 || eval_counter[] % trace_every == 0)
            pnow = metrics.p
            if !isempty(pnow)
                @printf("eval %4d  obj=%.8g  nχ²=%.8g  mχ²=%.8g  gzz=%.5g J1=%.5g J2=%.5g sigg=%.5g sigJ1=%.5g sigJ2=%.5g gzz2=%.5g sigg2=%.5g gperp=%.5g gperp2=%.5g r2=%.5g\n",
                    eval_counter[], metrics.total, metrics.neutron.redchi2, metrics.magnetization.redchi2,
                    pnow[:gzz], pnow[:J1_meV], pnow[:J2_meV], pnow[:sigma_gzz], pnow[:sigma_J1], pnow[:sigma_J2],
                    pnow[:gzz2], pnow[:sigma_gzz2], get(pnow, :gperp, NaN), get(pnow, :gperp2, NaN), shared_r2_from_cofit(pnow))
            else
                @printf("eval %4d  obj=%.8g\n", eval_counter[], metrics.total)
            end
        end
        return metrics.total
    end

    optres = nothing
    pbest = NF.unpack_unconstrained(u0, specs2)
    if run_optimization
        println("Starting Optim.jl Nelder-Mead co-fit in bounded-transformed parameters...")
        optres = optimize(obj, u0, NelderMead(), Optim.Options(iterations=maxiters, show_trace=show_trace))
        ubest = Optim.minimizer(optres)
        pbest = NF.unpack_unconstrained(ubest, specs2)
    else
        println("Skipping Optim.jl refinement; using the initial parameter values as the fixed model parameters.")
        println("Set run_optimization=true to refine these parameters again.")
    end

    pn_best = neutron_param_dict_from_cofit(pbest)
    r2_shared_final = shared_r2_from_cofit(pbest)

    println()
    if optres !== nothing
        println("Best optimizer result:")
        println(optres)
        println()
    end
    println(run_optimization ? "Best co-fit parameters:" : "Fixed model parameters:")
    for nm in [:gzz, :J1_meV, :J2_meV, :sigma_gzz, :sigma_J1, :sigma_J2,
               :gzz2, :sigma_gzz2, :gperp, :gperp2, :chi_vv_muB_per_T,
               :log10_second_kernel_relative_intensity]
        println(@sprintf("  %-56s % .10g", String(nm), pbest[nm]))
    end
    println(@sprintf("  %-56s % .10g", "r2_shared", r2_shared_final))

    println()
    println("Rerunning final neutron two-kernel model with final_n_samples_per_cut = ", final_n_samples_per_cut)
    neutron_final = evaluate_neutron_fit(pbest, data_scans;
        fields_T=fields_T,
        qtags=qtags,
        fit_windows_by_q=fit_windows_by_q,
        n_samples_per_cut=final_n_samples_per_cut,
        seed=seed,
        use_errors=use_errors,
        error_floor=error_floor,
        lattice=lattice,
        resolution=resolution,
        cuts=cuts,
        S=S,
        J_units=J_units,
        correlate_J1_J2=correlate_J1_J2,
        use_form_factor=use_form_factor,
        include_j2_formfactor=include_j2_formfactor,
        include_kfki=include_kfki,
        polarization=polarization,
    )

    println("Rerunning final magnetization model with final_n_samples_magnetization = ", final_n_samples_magnetization)
    mag_draws_final = make_standard_normal_draws_cofit(final_n_samples_magnetization; seed=mag_seed + 1)
    mag_final = evaluate_magnetization_fit(pbest, mag_fit, mag_draws_final, mag_quad;
        S_disp=S,
        S_imp=S_imp,
        temperature_K=temperature_K,
        dmax_mode=dmax_mode,
        response_mode=response_mode,
        magnetization_sigma_muB=magnetization_sigma_muB,
        magnetization_error_fraction_of_range=magnetization_error_fraction_of_range,
        magnetization_error_floor_muB=magnetization_error_floor_muB,
    )

    Bplot = collect(Float64(Bmin_plot_T):Float64(dB_plot_T):Float64(Bmax_plot_T))
    mag_curve = scaled_magnetization_curve(pbest, Bplot, mag_final.scale, mag_draws_final, mag_quad;
        S_disp=S,
        S_imp=S_imp,
        temperature_K=temperature_K,
        dmax_mode=dmax_mode,
        response_mode=response_mode,
    )

    comparison_data_scans = nothing
    comparison_model = nothing
    if add_14T_4p65_comparison_row
        println()
        println(@sprintf("Loading comparison neutron scans for Ei=%.5g meV, T=%.5g K", Float64(comparison_Ei_meV), Float64(comparison_temperature_K)))
        comparison_data_scans, _, _ = load_neutron_fit_data_1d_filtered(;
            base_dir=base_dir,
            data_mode=data_mode,
            Ei_meV=comparison_Ei_meV,
            temperature_K=comparison_temperature_K,
            background_kwargs...
        )
        println("Rerunning comparison neutron model for fields: ", collect(comparison_fields_T))
        comparison_model = evaluate_neutron_fit(pbest, comparison_data_scans;
            fields_T=comparison_fields_T,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            n_samples_per_cut=final_n_samples_per_cut,
            seed=seed + 14,
            use_errors=use_errors,
            error_floor=error_floor,
            lattice=lattice,
            resolution=resolution,
            cuts=cuts,
            S=S,
            J_units=J_units,
            correlate_J1_J2=correlate_J1_J2,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
        )
    end

    println()
    println(@sprintf("Final neutron reduced chi^2       = %.8g", neutron_final.redchi2))
    println(@sprintf("Final magnetization reduced chi^2 = %.8g", mag_final.redchi2))
    println(@sprintf("Final weighted objective          = %.8g",
        Float64(neutron_weight) * neutron_final.redchi2 + Float64(magnetization_weight) * mag_final.redchi2))
    println(@sprintf("Final neutron scale               = %.8g", neutron_final.scale))
    println(@sprintf("Final magnetization scale         = %.8g", mag_final.scale))
    println()

    summary_path = joinpath(outdir, "YZGO_neutron_magnetization_shared_fraction_cofit_summary.txt")
    neutron_fit_points_path = joinpath(outdir, "YZGO_neutron_magnetization_shared_fraction_cofit_neutron_fit_points.csv")
    neutron_model_csv_path = joinpath(outdir, "YZGO_neutron_magnetization_shared_fraction_cofit_neutron_scaled_model.csv")
    mag_model_csv_path = joinpath(outdir, "YZGO_neutron_magnetization_shared_fraction_cofit_magnetization_model.csv")
    combo_plot_path = joinpath(outdir, "YZGO_neutron_magnetization_shared_fraction_cofit_combo.png")

    write_cofit_summary(summary_path;
        specs=specs2,
        params=pbest,
        neutron_weight=neutron_weight,
        magnetization_weight=magnetization_weight,
        neutron_result=neutron_final,
        magnetization_result=mag_final,
        fields_T=fields_T,
        qtags=qtags,
        fit_windows_by_q=fit_windows_by_q,
        Bmin_fit_T=Bmin_fit_T,
        Bmax_fit_T=Bmax_fit_T,
        magnetization_sigma_muB=magnetization_sigma_muB,
        magnetization_error_fraction_of_range=magnetization_error_fraction_of_range,
        n_samples_per_cut=n_samples_per_cut,
        final_n_samples_per_cut=final_n_samples_per_cut,
        n_samples_magnetization=n_samples_magnetization,
        seed=seed,
        mag_seed=mag_seed,
        optimizer_result=optres,
    )
    NF.write_fit_points_csv(neutron_fit_points_path, neutron_final.bundle, neutron_final.scale, neutron_final.r2)
    NF.write_scaled_model_csv(neutron_model_csv_path, neutron_final.model;
        scale=neutron_final.scale,
        r2=neutron_final.r2,
        qtags=qtags,
    )
    write_magnetization_model_csv(mag_model_csv_path, mag_curve)

    fig = nothing
    if make_plots
        fig = plot_yzgo_cofit_combo(data_scans, neutron_final.model, mag_all, mag_fit, mag_curve;
            neutron_scale=neutron_final.scale,
            neutron_r2=neutron_final.r2,
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=fit_windows_by_q,
            neutron_ylims=neutron_plot_ylims,
            extra_data_scans=add_14T_4p65_comparison_row ? comparison_data_scans : nothing,
            extra_model_results_by_field=(add_14T_4p65_comparison_row && comparison_model !== nothing) ? comparison_model.model : nothing,
            extra_fields_T=add_14T_4p65_comparison_row ? comparison_fields_T : Float64[],
            extra_neutron_scale=neutron_final.scale,
            extra_neutron_r2=neutron_final.r2,
            extra_row_label=@sprintf("Ei=%.2f meV", Float64(comparison_Ei_meV)),
            outpath=combo_plot_path,
            save_png=true,
            display_fig=display_figures,
        )
    end

    println("Wrote:")
    println("  ", summary_path)
    println("  ", neutron_fit_points_path)
    println("  ", neutron_model_csv_path)
    println("  ", mag_model_csv_path)
    if make_plots
        println("  ", combo_plot_path)
    end

    return (; optimizer_result=optres,
              params=pbest,
              neutron=neutron_final,
              magnetization=mag_final,
              magnetization_curve=mag_curve,
              data_scans=data_scans,
              scans_raw=scans_raw,
              background_models=background_models,
              magnetization_data_all=mag_all,
              magnetization_data_fit=mag_fit,
              figure=fig,
              comparison_neutron=comparison_model,
              specs=specs2,
              paths=(summary=summary_path,
                     neutron_fit_points=neutron_fit_points_path,
                     neutron_model=neutron_model_csv_path,
                     magnetization_model=mag_model_csv_path,
                     combo_plot=combo_plot_path),
              weights=(neutron=neutron_weight, magnetization=magnetization_weight))
end



# =============================================================================
# No-optimization 2D comparison plot for fixed latest co-fit parameters
# =============================================================================
#
# This driver uses the model/simulation machinery above and the 2D data-loading
# style from plot_yzgo_2d_combo.jl.  It does NOT run any optimization.  It loads
# the Ei=4.65 meV, T=0.07 K, 9 T and 14 T 2D scans, simulates the preliminary
# best-fit two-component model, and makes one 2 x 2 plot per leg:
#
#   top row    = experimental data
#   bottom row = model from the fixed preliminary best-fit parameters
#   columns    = 9 T and 14 T
#
# To run:
#   julia plot_yzgo_2d_data_vs_prelim_model.jl
#
# Optional environment variables:
#   YZGO_DATA_DIR     folder containing yzgo_4p65meV_0p07K_9T_2d_leg*_SYM.dat
#   YZGO_OUT_DIR      output folder
#   MAKIE_BACKEND     GLMakie or CairoMakie
#
# Notes:
#   * The preliminary fit did not provide an absolute 2D neutron scale.  By
#     default, this script uses one analytic least-squares scale shared across
#     all plotted 2D panels.  Set MODEL_SCALE_MODE = :none or :manual if desired.
#   * The transverse g factors enter as neutron intensity prefactors gperp^2.
#   * The plotted model is fixed; no parameters are refined here.

# -----------------------------
# Plotting backend and controls
# -----------------------------

const YZGO_2D_DATA_DIR = get(ENV, "YZGO_DATA_DIR", raw"C:\Users\vdp\ORNL Dropbox\Daniel Pajerowski\YZGO\CNCS_data\SYM_2d_scans")
const YZGO_2D_OUT_DIR  = get(ENV, "YZGO_OUT_DIR", joinpath(YZGO_2D_DATA_DIR, "plots_2d_prelim_model_compare"))

const YZGO_2D_EI_TAG = "4p65"
const YZGO_2D_EI_MEV = 4.65
const YZGO_2D_TEMP_TAG = "0p07K"
const YZGO_2D_FIELDS_T = [9, 14]
const YZGO_2D_LEGS = [1]  # leg 1 only; leg 2 is symmetry-equivalent for this diagnostic plot

# :global_least_squares -> one scale shared across all plotted fields/legs
# :panel_least_squares  -> independent data/model scale for each panel
# :manual               -> use MANUAL_MODEL_SCALE
# :none                 -> no additional scaling
const MODEL_SCALE_MODE = :global_least_squares
const MANUAL_MODEL_SCALE = 0.04480834274  # latest 1D co-fit neutron_global_scale

# Simulation sampling. Increase for smoother final model maps.
const SAMPLES_PER_U_BIN_2D = 10_000
const MODEL_SEED_2D = 20260604

# Display controls.
const COLOR_MODE_2D = :shared_by_kind  # :shared_by_kind, :global, or :per_panel

# Experimental color scaling: the elastic line / near-elastic tail can dominate
# the heatmap contrast.  The data color range is therefore estimated from a
# restricted positive-energy window and a lower high-quantile than the model.
# This only changes plotting, not the model or any fit.
const DATA_COLOR_ENERGY_MIN_2D = 0.25
const DATA_COLOR_ENERGY_MAX_2D = Inf
const DATA_CLIP_HIGH_QUANTILE_2D = 0.950
const MODEL_CLIP_HIGH_QUANTILE_2D = 0.995
const CLIP_HIGH_QUANTILE_2D = 0.995  # fallback for generic arrays
const CLIP_LOW_QUANTILE_2D = 0.01
const FORCE_LO_ZERO_2D = true
const MASK_ZERO_INTENSITY_2D = false
const COLORMAP_2D = :viridis
const FIGURE_SIZE_2D = (1180, 850)

# Path construction.  :manual uses the same explicit leg definitions as the
# original 2D plotting/simulation convention; :header attempts to parse the
# exported .dat header.  If the simulated dispersion looks shifted along the
# path, toggling this is the first thing to check.
const PATH_SPEC_MODE_2D = :old_default_exact  # use the exact old analytic demo 2D path
const VERBOSE_PATH_SPECS_2D = true
const SAVE_PNG_2D = true
const SAVE_PDF_2D = lowercase(get(ENV, "MAKIE_BACKEND", "GLMakie")) == "cairomakie"
const SHOW_VERTICAL_GUIDES_2D = true
const GUIDE_XS_2D = [-1/3, 0.0, 1/3, 2/3, 1.0]
const ENERGY_LIMS_2D = (0.20, 3.20)  # avoid elastic-line-dominated display by default
const X_LIMS_2D = nothing       # e.g. (-1/3, 1.0)

# Latest co-fit best parameters used as the fixed 2D model.
# gperp_ratio is a fitted neutron-only relative transverse intensity prefactor.
# It is not an absolute physical g_perp.  Because the model also has an overall
# neutron intensity scale, the absolute dispersive transverse prefactor is a
# gauge convention and only the flat/dispersive ratio is identifiable.

# ---------------------------------------------------------------------------
# Latest co-fit parameters for 2D model comparison
# ---------------------------------------------------------------------------
#
# Repo refactor note:
#
# Older versions of this script carried a hard-coded NamedTuple of the latest
# co-fit parameters. The repo version instead reads:
#
#   configs/best_fit_parameters.toml
#
# via ENV["YZGO_BEST_FIT_PARAMETERS_TOML"], which is set by
# scripts/plot_2d_data_vs_model.jl before this file is included.
#
# This makes the 1D co-fit and the 2D visualization use the same parameter
# source.

using TOML

function _toml_lookup_2d(config::Dict, sections::Vector{String}, key::String)
    for section in sections
        if haskey(config, section)
            table = config[section]
            if table isa Dict && haskey(table, key)
                return table[key]
            end
        end
    end

    tried = join(["[$s].$key" for s in sections], ", ")
    error("Missing required TOML value for 2D model parameters. Tried: $tried")
end

function _toml_positive_or_log10_2d(
    config::Dict;
    sections::Vector{String},
    value_key::String,
    log10_key::String,
)
    for section in sections
        if !haskey(config, section)
            continue
        end

        table = config[section]
        if !(table isa Dict)
            continue
        end

        if haskey(table, value_key)
            return Float64(table[value_key])
        end

        if haskey(table, log10_key)
            return 10.0 ^ Float64(table[log10_key])
        end
    end

    tried_value = join(["[$s].$value_key" for s in sections], ", ")
    tried_log10 = join(["[$s].$log10_key" for s in sections], ", ")
    error("Missing required TOML value for 2D model parameters. Tried: $tried_value or $tried_log10")
end

function load_2d_model_parameters_from_best_fit_toml(path::AbstractString)
    if !isfile(path)
        error("Could not find best-fit parameter TOML for 2D model: $path")
    end

    config = TOML.parsefile(path)

    physical_sections = ["physical", "initial_guess", "parameters"]

    extrinsic_sections = [
        "neutron_extrinsic",
        "extrinsic",
        "physical",
        "initial_guess",
        "parameters",
    ]

    sigma_J = Float64(_toml_lookup_2d(config, physical_sections, "sigma_J"))

    r2_shared = _toml_positive_or_log10_2d(
        config;
        sections = extrinsic_sections,
        value_key = "second_kernel_relative_intensity",
        log10_key = "log10_second_kernel_relative_intensity",
    )

    neutron_global_scale = _toml_positive_or_log10_2d(
        config;
        sections = extrinsic_sections,
        value_key = "neutron_global_scale",
        log10_key = "log10_neutron_scale",
    )

    gperp_ratio = Float64(_toml_lookup_2d(config, physical_sections, "gperp_ratio"))

    return (;
        gzz = Float64(_toml_lookup_2d(config, physical_sections, "gzz")),
        J1_meV = Float64(_toml_lookup_2d(config, physical_sections, "J1_meV")),
        J2_meV = Float64(_toml_lookup_2d(config, physical_sections, "J2_meV")),

        # Latest model convention: one shared fractional exchange-disorder width.
        sigma_gzz = Float64(_toml_lookup_2d(config, physical_sections, "sigma_gzz")),
        sigma_J = sigma_J,

        gzz2 = Float64(_toml_lookup_2d(config, physical_sections, "gzz2")),
        sigma_gzz2 = Float64(_toml_lookup_2d(config, physical_sections, "sigma_gzz2")),

        # Latest model convention: one effective neutron transverse-intensity ratio.
        # This is the relative matrix-element/intensity amplitude of the flat
        # component compared with the dispersive component, not a physical g⊥.
        gperp_ratio = gperp_ratio,

        chi_vv_muB_per_T = Float64(_toml_lookup_2d(config, physical_sections, "chi_vv_muB_per_T")),

        # Positive extrinsic neutron intensity parameters.
        second_kernel_relative_intensity = r2_shared,
        neutron_global_scale = neutron_global_scale,
    )
end

const BEST_FIT_PARAMETERS_TOML_FOR_2D = get(
    ENV,
    "YZGO_BEST_FIT_PARAMETERS_TOML",
    "",
)

if isempty(BEST_FIT_PARAMETERS_TOML_FOR_2D)
    error(
        "ENV[\"YZGO_BEST_FIT_PARAMETERS_TOML\"] is not set. " *
        "Run this through scripts/plot_2d_data_vs_model.jl."
    )
end

const LAST_COFIT_BEST_PARAMS_2D = load_2d_model_parameters_from_best_fit_toml(
    BEST_FIT_PARAMETERS_TOML_FOR_2D,
)

# Canonical parameter dictionary used by the 2D plotting layer.
#
# User-facing/science-facing keys should stay in the latest model language:
#
#   sigma_J       one shared fractional exchange-disorder width
#   gperp_ratio   relative transverse neutron matrix-element amplitude
#
# Older helper functions in the remaining legacy analytical engine still expect
# :sigma_J1/:sigma_J2 and :gperp/:gperp2.  Do not expose those as model
# parameters; use `_nf_legacy_fit_dict_for_2d` only at the boundary where we
# call those helper functions.
const PRELIM_BEST_PARAMS = Dict{Symbol,Float64}(
    :gzz => LAST_COFIT_BEST_PARAMS_2D.gzz,
    :J1_meV => LAST_COFIT_BEST_PARAMS_2D.J1_meV,
    :J2_meV => LAST_COFIT_BEST_PARAMS_2D.J2_meV,
    :sigma_gzz => LAST_COFIT_BEST_PARAMS_2D.sigma_gzz,
    :sigma_J => LAST_COFIT_BEST_PARAMS_2D.sigma_J,
    :gzz2 => LAST_COFIT_BEST_PARAMS_2D.gzz2,
    :sigma_gzz2 => LAST_COFIT_BEST_PARAMS_2D.sigma_gzz2,
    :gperp_ratio => LAST_COFIT_BEST_PARAMS_2D.gperp_ratio,
    :chi_vv_muB_per_T => LAST_COFIT_BEST_PARAMS_2D.chi_vv_muB_per_T,
    :second_kernel_relative_intensity => LAST_COFIT_BEST_PARAMS_2D.second_kernel_relative_intensity,
    :neutron_global_scale => LAST_COFIT_BEST_PARAMS_2D.neutron_global_scale,
)

function _nf_legacy_fit_dict_for_2d(p::Dict{Symbol,Float64})
    # Boundary adapter only. This is not a second parameterization.
    #
    # The analytical helper functions in NF still use older dictionary keys.
    # The current model has one sigma_J and one gperp_ratio, so this adapter
    # maps:
    #
    #   sigma_J -> sigma_J1 = sigma_J2
    #   gperp_ratio -> gperp = 1, gperp2 = gperp_ratio
    #
    return Dict{Symbol,Float64}(
        :gzz => p[:gzz],
        :J1_meV => p[:J1_meV],
        :J2_meV => p[:J2_meV],
        :sigma_gzz => p[:sigma_gzz],
        :sigma_J1 => p[:sigma_J],
        :sigma_J2 => p[:sigma_J],
        :gzz2 => p[:gzz2],
        :sigma_gzz2 => p[:sigma_gzz2],
        :gperp => 1.0,
        :gperp2 => p[:gperp_ratio],
        :log10_second_kernel_relative_intensity => log10(p[:second_kernel_relative_intensity]),
        :log10_scale => log10(p[:neutron_global_scale]),
        :neutron_global_scale => p[:neutron_global_scale],
        :chi_vv_muB_per_T => p[:chi_vv_muB_per_T],
    )
end

function _neutron_intensity_scales_from_canonical_2d(p::Dict{Symbol,Float64})
    return (;
        dispersive = 1.0,
        nondispersive = max(p[:gperp_ratio], 0.0)^2,
    )
end

# -----------------------------
# 2D data structure/helpers
# -----------------------------

struct Scan2DCompare
    file::String
    header::String
    xlabel::String
    x::Vector{Float64}
    e::Vector{Float64}
    z::Matrix{Float64}
end

function _header_line_2d(file::AbstractString)
    open(file, "r") do io
        return strip(readline(io))
    end
end

function _path_label_from_header_2d(header::AbstractString)
    m = match(r"Error\s+(.+?)\s+DeltaE", header)
    return m === nothing ? "path coordinate (rlu)" : "path coordinate " * strip(m.captures[1]) * " (rlu)"
end

function _find_scan_file_2d(data_dir::AbstractString, field_T::Integer, leg::Integer)
    pattern = Regex("^yzgo_$(YZGO_2D_EI_TAG)meV_$(YZGO_2D_TEMP_TAG)_$(field_T)T_2d_leg$(leg)_SYM\\.dat\$")
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

function _read_scan2d_compare(file::AbstractString; mask_zero::Bool=false)
    header = _header_line_2d(file)
    data = DelimitedFiles.readdlm(file, Float64; comments=true, comment_char='#')
    size(data, 2) >= 4 || error("Expected at least 4 numeric columns in $(file), got $(size(data,2)).")

    intensity = data[:, 1]
    xcol = data[:, 3]
    ecol = data[:, 4]
    xs = sort(collect(unique(xcol)))
    es = sort(collect(unique(ecol)))
    xindex = Dict(v => i for (i, v) in enumerate(xs))
    eindex = Dict(v => i for (i, v) in enumerate(es))
    z = fill(NaN, length(xs), length(es))
    for r in axes(data, 1)
        val = intensity[r]
        if mask_zero && iszero(val)
            val = NaN
        end
        z[xindex[xcol[r]], eindex[ecol[r]]] = val
    end
    return Scan2DCompare(file, header, _path_label_from_header_2d(header), xs, es, z)
end

function _finite_values_2d(z)
    vals = Float64[]
    for v in z
        if isfinite(v)
            push!(vals, Float64(v))
        end
    end
    return vals
end

function _robust_colorrange_from_values_2d(vals::Vector{Float64}; high_quantile::Real=CLIP_HIGH_QUANTILE_2D)
    if isempty(vals)
        return (0.0, 1.0)
    end
    lo = FORCE_LO_ZERO_2D ? 0.0 : quantile(vals, CLIP_LOW_QUANTILE_2D)
    hi = quantile(vals, Float64(high_quantile))
    if !isfinite(lo) || !isfinite(hi) || hi <= lo
        lo = FORCE_LO_ZERO_2D ? 0.0 : minimum(vals)
        hi = maximum(vals)
    end
    if hi <= lo
        hi = lo + 1.0
    end
    return (lo, hi)
end

function _robust_colorrange_2d(arrays::Vector{Matrix{Float64}}; high_quantile::Real=CLIP_HIGH_QUANTILE_2D)
    vals = Float64[]
    for z in arrays
        append!(vals, _finite_values_2d(z))
    end
    return _robust_colorrange_from_values_2d(vals; high_quantile=high_quantile)
end

function _finite_values_energy_window_2d(scan::Scan2DCompare; emin::Real=DATA_COLOR_ENERGY_MIN_2D, emax::Real=DATA_COLOR_ENERGY_MAX_2D)
    vals = Float64[]
    for (ie, E) in enumerate(scan.e)
        if Float64(emin) <= E <= Float64(emax)
            for ix in eachindex(scan.x)
                v = scan.z[ix, ie]
                if isfinite(v)
                    push!(vals, Float64(v))
                end
            end
        end
    end
    return vals
end

function _data_colorrange_2d(scans::Vector{Scan2DCompare})
    vals = Float64[]
    for s in scans
        append!(vals, _finite_values_energy_window_2d(s))
    end
    if isempty(vals)
        @warn "No finite values in requested data color energy window; falling back to all energies" emin=DATA_COLOR_ENERGY_MIN_2D emax=DATA_COLOR_ENERGY_MAX_2D
        for s in scans
            append!(vals, _finite_values_2d(s.z))
        end
    end
    return _robust_colorrange_from_values_2d(vals; high_quantile=DATA_CLIP_HIGH_QUANTILE_2D)
end

function _step_from_centers(v::Vector{Float64}; fallback::Float64)
    length(v) >= 2 || return fallback
    ds = diff(v)
    return median(ds)
end

# -----------------------------
# Header-to-model path helpers
# -----------------------------

function _split_bracket_vectors(header::AbstractString)
    # Returns strings inside [ ... ]. For the expected header this is:
    #   path expression, integration expression, integration expression
    return [m.captures[1] for m in eachmatch(r"\[([^\]]+)\]", header)]
end

function _parse_component_coeff(comp::AbstractString)
    s = replace(strip(comp), " " => "")
    isempty(s) && return 0.0
    s == "0" && return 0.0
    s == "H" && return 1.0
    s == "K" && return 1.0
    s == "L" && return 1.0
    s == "-H" && return -1.0
    s == "-K" && return -1.0
    s == "-L" && return -1.0

    # Handles strings like -0.5H, 0.5K, 2L. Constant offsets are intentionally
    # not supported by PathCutSpec2D; warn if one is encountered elsewhere.
    m = match(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))\*?([HKL])$", s)
    if m !== nothing
        return parse(Float64, m.captures[1])
    end
    try
        return parse(Float64, s)
    catch
        @warn "Could not parse vector component in 2D header; treating as 0" comp
        return 0.0
    end
end

function _parse_vector_expr(expr::AbstractString)
    parts = split(expr, ",")
    if length(parts) != 3
        @warn "Expected a 3-component bracket vector; using zeros" expr
        return (0.0, 0.0, 0.0)
    end
    return tuple((_parse_component_coeff(p) for p in parts)...)
end

function _default_range_for_basis(vec::NTuple{3,Float64})
    # Heuristics matching the previously used 2D CNCS cuts:
    #   [0,0,L] integrated roughly around L=0
    #   [0,K,0] integrated around K=0.5 for the leg-1 path
    #   other in-plane perpendicular integrations around 0
    if abs(vec[1]) < 1e-12 && abs(vec[2]) < 1e-12 && abs(vec[3]) > 1e-12
        return (-0.25, 0.25)
    elseif abs(vec[1]) < 1e-12 && abs(vec[2] - 1.0) < 1e-12 && abs(vec[3]) < 1e-12
        return (0.45, 0.55)
    else
        return (-0.05, 0.05)
    end
end


function _old_default_exact_spec_2d()
    # This is intentionally the exact PathCutSpec2D convention from
    # YZGO_analytic_1d2d_glmakie_stratified_1d2d.jl:default_leg_cut_2d().
    # It does NOT infer u_range or binning from the exported 2D data file.
    # This avoids a subtle but important mismatch in which the newer script
    # tried to reinterpret the experimental grid/header as the model path.
    return NF.PathCutSpec2D(
        name = "old_default_K1_to_Gamma1_leg",
        uvec = (1.0, -0.5, 0.0),
        vvec = (0.0, 1.0, 0.0),
        wvec = (0.0, 0.0, 1.0),
        u_range = (-1.0 / 3.0, 1.0),
        v_range = (0.45, 0.55),
        w_range = (-0.25, 0.25),
        du = 0.01,
        E_range = (0.0, 4.0),
        dE_meV = 0.05,
    )
end

function _old_analytic_path_for_leg_2d(scan::Scan2DCompare, leg::Integer)
    # Reproduce the older analytic 2D simulation convention as closely as
    # possible.  This treats the x coordinate in the exported 2D cut as the
    # path coordinate u and explicitly defines the reciprocal-space path and
    # integration directions, rather than trying to infer them from the header.
    #
    # Leg 1 follows the old demo/default path:
    #   Q(u, v, w) = u * [1, -1/2, 0] + v * [0, 1, 0] + w * [0, 0, 1]
    # with v integrated around 1/2.  This gives Q = (u, 1/2 - u/2, L), i.e.
    # K1 -> M1 -> K -> Gamma1 when u goes from -1/3 to 1.
    #
    # Leg 2 uses the corresponding header convention:
    #   Q(u, v, w) = u * [1, 1, 0] + v * [1, -1, 0] + w * [0, 0, 1]
    # with v integrated around 0.
    if leg == 1
        uvec = (1.0, -0.5, 0.0)
        vvec = (0.0, 1.0, 0.0)
        wvec = (0.0, 0.0, 1.0)
        v_range = (0.45, 0.55)
        w_range = (-0.25, 0.25)
    elseif leg == 2
        uvec = (1.0, 1.0, 0.0)
        vvec = (1.0, -1.0, 0.0)
        wvec = (0.0, 0.0, 1.0)
        v_range = (-0.05, 0.05)
        w_range = (-0.25, 0.25)
    else
        error("No old analytic 2D path convention defined for leg=$(leg)")
    end
    return (; uvec, vvec, wvec, v_range, w_range)
end

function _manual_vectors_for_leg_2d(leg::Integer)
    if leg == 1
        # Header convention: [H,-0.5H,0] DeltaE [0,0,L] [0,K,0].
        # With K integrated near 0.5 this traces (H, 0.5 - H/2, L).
        return ((1.0, -0.5, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    elseif leg == 2
        # Header convention: [H,H,0] DeltaE [0,0,L] [K,-K,0].
        return ((1.0, 1.0, 0.0), (1.0, -1.0, 0.0), (0.0, 0.0, 1.0))
    else
        error("No manual 2D path convention defined for leg=$(leg)")
    end
end

function _spec_from_scan(scan::Scan2DCompare; leg::Integer)
    if PATH_SPEC_MODE_2D == :old_default_exact
        leg == 1 || error("PATH_SPEC_MODE_2D=:old_default_exact is defined for leg 1 only; leg 2 is symmetry-equivalent and disabled by default.")
        spec = _old_default_exact_spec_2d()
        if VERBOSE_PATH_SPECS_2D
            @info "2D model path spec" leg mode=PATH_SPEC_MODE_2D header=scan.header uvec=spec.uvec vvec=spec.vvec wvec=spec.wvec u_range=spec.u_range v_range=spec.v_range w_range=spec.w_range E_range=spec.E_range du=spec.du dE_meV=spec.dE_meV
        end
        return spec
    end

    v_range_override = nothing
    w_range_override = nothing
    if PATH_SPEC_MODE_2D == :old_analytic
        oldpath = _old_analytic_path_for_leg_2d(scan, leg)
        uvec = oldpath.uvec
        vvec = oldpath.vvec
        wvec = oldpath.wvec
        v_range_override = oldpath.v_range
        w_range_override = oldpath.w_range
    elseif PATH_SPEC_MODE_2D == :manual
        uvec, vvec, wvec = _manual_vectors_for_leg_2d(leg)
    elseif PATH_SPEC_MODE_2D == :header
        brackets = _split_bracket_vectors(scan.header)
        if length(brackets) >= 3
            uvec = _parse_vector_expr(brackets[1])
            vvec = _parse_vector_expr(brackets[3])  # keep [0,K,0] as v when present
            wvec = _parse_vector_expr(brackets[2])  # keep [0,0,L] as w when present
        else
            @warn "Could not parse path/integration vectors from header; using old analytic fallback for leg" leg header=scan.header
            oldpath = _old_analytic_path_for_leg_2d(scan, leg)
            uvec = oldpath.uvec
            vvec = oldpath.vvec
            wvec = oldpath.wvec
            v_range_override = oldpath.v_range
            w_range_override = oldpath.w_range
        end
    else
        error("Unknown PATH_SPEC_MODE_2D = $(PATH_SPEC_MODE_2D)")
    end

    du = _step_from_centers(scan.x; fallback=0.01)
    dE = _step_from_centers(scan.e; fallback=0.05)
    u_range = (minimum(scan.x) - du/2, maximum(scan.x) + du/2)
    E_range = (minimum(scan.e) - dE/2, maximum(scan.e) + dE/2)
    v_range = v_range_override === nothing ? _default_range_for_basis(vvec) : v_range_override
    w_range = w_range_override === nothing ? _default_range_for_basis(wvec) : w_range_override

    spec = NF.PathCutSpec2D(
        name = "leg$(leg)_$(PATH_SPEC_MODE_2D)",
        uvec = uvec,
        vvec = vvec,
        wvec = wvec,
        u_range = u_range,
        v_range = v_range,
        w_range = w_range,
        du = du,
        E_range = E_range,
        dE_meV = dE,
    )
    if VERBOSE_PATH_SPECS_2D
        @info "2D model path spec" leg mode=PATH_SPEC_MODE_2D header=scan.header uvec=spec.uvec vvec=spec.vvec wvec=spec.wvec u_range=spec.u_range v_range=spec.v_range w_range=spec.w_range E_range=spec.E_range du=spec.du dE_meV=spec.dE_meV
    end
    return spec
end

# -----------------------------
# Fixed preliminary model simulation
# -----------------------------

function _model_params_for_field_2d(p::Dict{Symbol,Float64}, field_T::Real)
    p_nf = _nf_legacy_fit_dict_for_2d(p)

    disp = NF._dispersive_model_params_from_fit_dict(p_nf, field_T; S=0.5)
    flat = NF._nondispersive_model_params_from_fit_dict(p_nf, field_T; S=0.5)
    ddisp = NF._dispersive_disorder_params_from_fit_dict(p_nf; J_units=:fractional, correlate_J1_J2=false)
    dflat = NF._nondispersive_disorder_params_from_fit_dict(p_nf)
    gscales = _neutron_intensity_scales_from_canonical_2d(p)

    return (; disp, flat, ddisp, dflat, gscales, r2=p[:second_kernel_relative_intensity])
end

function _simulate_prelim_model_2d(scan::Scan2DCompare, field_T::Real, leg::Integer;
                                   samples_per_u_bin::Int=SAMPLES_PER_U_BIN_2D,
                                   seed::Int=MODEL_SEED_2D)
    p = PRELIM_BEST_PARAMS
    mp = _model_params_for_field_2d(p, field_T)
    spec = _spec_from_scan(scan; leg=leg)
    lattice = NF.demo_defaults().lattice
    resolution = NF.ResolutionParams(Ei_meV=YZGO_2D_EI_MEV)
    field_seed = seed + 10_000 * leg + round(Int, 1000 * Float64(field_T))

    disp = NF.simulate_path_map_2d_stratified(
        spec, mp.disp;
        lattice=lattice,
        disorder=mp.ddisp,
        resolution=resolution,
        rng=MersenneTwister(field_seed),
        samples_per_u_bin=samples_per_u_bin,
        use_form_factor=true,
        include_j2_formfactor=true,
        polarization=:transverse_c,
        include_kfki=true,
        intensity_scale=mp.gscales.dispersive,
    )
    flat = NF.simulate_path_map_2d_stratified(
        spec, mp.flat;
        lattice=lattice,
        disorder=mp.dflat,
        resolution=resolution,
        rng=MersenneTwister(field_seed + 7919),
        samples_per_u_bin=samples_per_u_bin,
        use_form_factor=true,
        include_j2_formfactor=true,
        polarization=:transverse_c,
        include_kfki=true,
        intensity_scale=mp.gscales.nondispersive,
    )

    zmodel = disp.intensity .+ mp.r2 .* flat.intensity
    return (; spec, disp, flat, x=disp.u_centers, e=disp.E_centers_meV, z=zmodel, r2=mp.r2)
end

function _load_all_2d_data()
    scans = Dict{Tuple{Int,Int},Scan2DCompare}()
    for leg in YZGO_2D_LEGS, field_T in YZGO_2D_FIELDS_T
        file = _find_scan_file_2d(YZGO_2D_DATA_DIR, field_T, leg)
        if file === nothing
            @warn "Missing requested 2D scan" data_dir=YZGO_2D_DATA_DIR field_T leg
            continue
        end
        scans[(leg, field_T)] = _read_scan2d_compare(file; mask_zero=MASK_ZERO_INTENSITY_2D)
    end
    return scans
end

function _simulate_all_2d_models(scans)
    models = Dict{Tuple{Int,Int},Any}()
    for ((leg, field_T), scan) in scans
        @info "Simulating preliminary 2D model" leg field_T samples_per_u_bin=SAMPLES_PER_U_BIN_2D
        models[(leg, field_T)] = _simulate_prelim_model_2d(scan, field_T, leg)
    end
    return models
end

function _least_squares_scale(data_arrays, model_arrays)
    # Robustly estimate one multiplicative model scale.  The experimental 2D
    # files do not always have exactly the same grid dimensions as the model
    # arrays, even when the model was generated from the scan metadata.  Rather
    # than requiring identical linear indices, compare the overlapping matrix
    # region for each data/model pair.  This is only used for display scaling;
    # no optimization is performed in this script.
    num = 0.0
    den = 0.0
    for (d, m) in zip(data_arrays, model_arrays)
        nx = min(size(d, 1), size(m, 1))
        ne = min(size(d, 2), size(m, 2))
        if size(d) != size(m)
            @warn "Data/model grid size mismatch during display scaling; using overlapping matrix region" data_size=size(d) model_size=size(m) overlap=(nx, ne)
        end
        for ix in 1:nx, ie in 1:ne
            di = d[ix, ie]
            mi = m[ix, ie]
            if isfinite(di) && isfinite(mi)
                num += di * mi
                den += mi * mi
            end
        end
    end
    return den > 0 ? max(0.0, num / den) : 1.0
end

function _scaled_models(scans, models)
    scale_by_key = Dict{Tuple{Int,Int},Float64}()
    if MODEL_SCALE_MODE == :none
        for k in keys(models); scale_by_key[k] = 1.0; end
    elseif MODEL_SCALE_MODE == :manual
        for k in keys(models); scale_by_key[k] = MANUAL_MODEL_SCALE; end
    elseif MODEL_SCALE_MODE == :panel_least_squares
        for k in keys(models)
            scale_by_key[k] = _least_squares_scale([scans[k].z], [models[k].z])
        end
    elseif MODEL_SCALE_MODE == :global_least_squares
        common_keys = sort(collect(intersect(keys(scans), keys(models))))
        s = _least_squares_scale([scans[k].z for k in common_keys], [models[k].z for k in common_keys])
        for k in keys(models); scale_by_key[k] = s; end
    else
        error("Unknown MODEL_SCALE_MODE = $(MODEL_SCALE_MODE)")
    end

    scaled = Dict{Tuple{Int,Int},Matrix{Float64}}()
    for k in keys(models)
        scaled[k] = scale_by_key[k] .* models[k].z
    end
    return scaled, scale_by_key
end

# -----------------------------
# Plotting
# -----------------------------

function _axis_limits_2d!(ax, x, e)
    if X_LIMS_2D === nothing
        xlims!(ax, minimum(x), maximum(x))
    else
        xlims!(ax, X_LIMS_2D...)
    end
    if ENERGY_LIMS_2D === nothing
        ylims!(ax, minimum(e), maximum(e))
    else
        ylims!(ax, ENERGY_LIMS_2D...)
    end
end

function _add_guides_2d!(ax, x)
    SHOW_VERTICAL_GUIDES_2D || return
    xmin, xmax = X_LIMS_2D === nothing ? (minimum(x), maximum(x)) : X_LIMS_2D
    guides = [xx for xx in GUIDE_XS_2D if xmin <= xx <= xmax]
    isempty(guides) || vlines!(ax, guides; color=(:white, 0.45), linewidth=1)
end

function _colorranges_2d_for_leg(scans, scaled_models, leg::Integer)
    keys_leg = [(leg, f) for f in YZGO_2D_FIELDS_T if haskey(scans, (leg, f))]
    if COLOR_MODE_2D == :global
        # Keep data and model colorbars separate even in global mode, because
        # near-elastic experimental leakage can otherwise wash out the model.
        data_cr = _data_colorrange_2d(collect(values(scans)))
        model_cr = _robust_colorrange_2d(collect(values(scaled_models)); high_quantile=MODEL_CLIP_HIGH_QUANTILE_2D)
        return Dict(:data => data_cr, :model => model_cr)
    elseif COLOR_MODE_2D == :shared_by_kind
        data_cr = _data_colorrange_2d([scans[k] for k in keys_leg])
        model_cr = _robust_colorrange_2d([scaled_models[k] for k in keys_leg if haskey(scaled_models, k)]; high_quantile=MODEL_CLIP_HIGH_QUANTILE_2D)
        return Dict(:data => data_cr, :model => model_cr)
    elseif COLOR_MODE_2D == :per_panel
        return Dict(:data => nothing, :model => nothing)
    else
        error("Unknown COLOR_MODE_2D = $(COLOR_MODE_2D)")
    end
end

function make_2d_prelim_data_model_plot_for_leg(leg::Integer, scans, models, scaled_models, scale_by_key)
    fig = Figure(size=FIGURE_SIZE_2D, fontsize=16)
    title = @sprintf("YZGO 2D data vs preliminary fixed model, leg %d, Ei=%.2f meV, T=0.07 K", leg, YZGO_2D_EI_MEV)
    Label(fig[0, 1:3], title, fontsize=21, font=:bold, tellwidth=false)

    cranges = _colorranges_2d_for_leg(scans, scaled_models, leg)
    data_hm = nothing
    model_hm = nothing
    first_xlabel = "Path coordinate (rlu)"

    for (icol, field_T) in enumerate(YZGO_2D_FIELDS_T)
        key = (leg, field_T)
        if !haskey(scans, key)
            for irow in 1:2
                ax = Axis(fig[irow, icol], title=irow==1 ? @sprintf("%d T", field_T) : "")
                hidedecorations!(ax); hidespines!(ax)
                text!(ax, 0.5, 0.5; text="missing", align=(:center, :center), space=:relative)
            end
            continue
        end

        scan = scans[key]
        model = models[key]
        zmodel = scaled_models[key]
        first_xlabel = scan.xlabel

        axd = Axis(fig[1, icol],
            title=@sprintf("%d T", field_T),
            ylabel=icol == 1 ? "Data\nΔE (meV)" : "",
            xlabel="",
        )
        dcr = COLOR_MODE_2D == :per_panel ? _data_colorrange_2d([scan]) : cranges[:data]
        data_hm = heatmap!(axd, scan.x, scan.e, scan.z;
            colormap=COLORMAP_2D,
            colorrange=dcr,
            nan_color=:lightgray,
        )
        _axis_limits_2d!(axd, scan.x, scan.e)
        _add_guides_2d!(axd, scan.x)

        axm = Axis(fig[2, icol],
            ylabel=icol == 1 ? "Model\nΔE (meV)" : "",
            xlabel=scan.xlabel,
            title=@sprintf("scale = %.4g", scale_by_key[key]),
        )
        mcr = COLOR_MODE_2D == :per_panel ? _robust_colorrange_2d([zmodel]; high_quantile=MODEL_CLIP_HIGH_QUANTILE_2D) : cranges[:model]
        model_hm = heatmap!(axm, model.x, model.e, zmodel;
            colormap=COLORMAP_2D,
            colorrange=mcr,
            nan_color=:lightgray,
        )
        _axis_limits_2d!(axm, model.x, model.e)
        _add_guides_2d!(axm, model.x)
    end

    if data_hm !== nothing
        Colorbar(fig[1, 3], data_hm; label="Data intensity")
    end
    if model_hm !== nothing
        Colorbar(fig[2, 3], model_hm; label="Scaled model intensity")
    end

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)
    return fig
end

function print_prelim_parameter_summary_2d()
    println("Fixed latest co-fit parameters used for 2D comparison")
    println("----------------------------------------------------------")
    for nm in [:gzz, :J1_meV, :J2_meV, :sigma_gzz, :sigma_J,
               :gzz2, :sigma_gzz2, :gperp_ratio, :chi_vv_muB_per_T,
               :second_kernel_relative_intensity, :neutron_global_scale]
        println(@sprintf("%-42s %.10g", String(nm), PRELIM_BEST_PARAMS[nm]))
    end
    gsc = _neutron_intensity_scales_from_canonical_2d(PRELIM_BEST_PARAMS)
    println(@sprintf("%-42s %.10g", "dispersive_transverse_intensity_factor", gsc.dispersive))
    println(@sprintf("%-42s %.10g", "flat_transverse_intensity_factor", gsc.nondispersive))
    println(@sprintf("%-42s %.10g", "manual_model_scale", MANUAL_MODEL_SCALE))
end

function run_yzgo_2d_prelim_model_comparison()
    @info "YZGO fixed-parameter 2D data/model comparison" data_dir=YZGO_2D_DATA_DIR out_dir=YZGO_2D_OUT_DIR backend=get(ENV, "MAKIE_BACKEND", "GLMakie") scale_mode=MODEL_SCALE_MODE path_mode=PATH_SPEC_MODE_2D data_color_window=(DATA_COLOR_ENERGY_MIN_2D, DATA_COLOR_ENERGY_MAX_2D) data_clip=DATA_CLIP_HIGH_QUANTILE_2D
    print_prelim_parameter_summary_2d()
    mkpath(YZGO_2D_OUT_DIR)

    scans = _load_all_2d_data()
    isempty(scans) && error("No requested 2D data files were found in $(YZGO_2D_DATA_DIR).")
    models = _simulate_all_2d_models(scans)
    scaled_models, scale_by_key = _scaled_models(scans, models)

    figs = Figure[]
    for leg in YZGO_2D_LEGS
        any(haskey(scans, (leg, f)) for f in YZGO_2D_FIELDS_T) || continue
        fig = make_2d_prelim_data_model_plot_for_leg(leg, scans, models, scaled_models, scale_by_key)
        push!(figs, fig)

        if SAVE_PNG_2D
            path = joinpath(YZGO_2D_OUT_DIR, @sprintf("YZGO_2d_data_vs_prelim_model_Ei4p65_oldpath_leg%d_9T14T.png", leg))
            save(path, fig)
            @info "Saved" path
        end
        if SAVE_PDF_2D
            path = joinpath(YZGO_2D_OUT_DIR, @sprintf("YZGO_2d_data_vs_prelim_model_Ei4p65_oldpath_leg%d_9T14T.pdf", leg))
            if get(ENV, "MAKIE_BACKEND", "GLMakie") == "CairoMakie"
                try
                    CairoMakie.activate!()
                    save(path, fig)
                    @info "Saved" path
                catch err
                    @warn "PDF save failed; PNG was already written" exception=(err, catch_backtrace())
                end
            end
        end
    end

    if lowercase(get(ENV, "MAKIE_BACKEND", "GLMakie")) != "cairomakie"
        display.(figs)
    end

    return (; scans, models, scaled_models, scale_by_key, figures=figs)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result_2d_prelim = run_yzgo_2d_prelim_model_comparison()
end

