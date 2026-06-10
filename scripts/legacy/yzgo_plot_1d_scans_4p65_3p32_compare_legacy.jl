# yzgo_plot_1d_scans_4p65_3p32_compare.jl
#
# Compare background-subtracted YZGO 1D energy scans for Ei = 4.65 meV and
# Ei = 3.32 meV.
#
# First-pass background rules:
#   Ei = 4.65 meV:
#       Keep the same tail / minimum-field / structured-residual background style
#       from the starter script yzgo_plot_1d_scans.jl.
#
#   Ei = 3.32 meV:
#       For each Q scan, use the point-by-point minimum intensity over these four
#       configurations as the background:
#           0.07 K, 0 T
#           0.07 K, 9 T
#           0.07 K, 14 T
#           20 K,   0 T
#
# Notes:
#   * There is no gamma-point scan for Ei = 3.32 meV, so the comparison plot leaves
#     that panel marked as not measured.
#   * The parser accepts both 0p07K and 20K filename temperature tokens.
#   * This script is standalone and does not require LsqFit/Interpolations.
#
# Usage option 1, from a Julia REPL:
#   include("yzgo_plot_1d_scans_4p65_3p32_compare.jl")
#   scans, bgsub_scans, bg_models, fig = main()
#
# Usage option 2, from a terminal:
#   julia yzgo_plot_1d_scans_4p65_3p32_compare.jl
#
# If your data are not in BASEDIR below, either edit BASEDIR or set:
#   ENV["YZGO_1D_SCAN_DIR"] = raw"C:\path\to\SYM_1d_scans"
# before calling main().
#
# Repo migration note:
#   This legacy implementation is intended to be called by
#   scripts/compare_1d_4p65_3p32_backgrounds.jl, which supplies repo-relative
#   paths and plotting controls. This file no longer runs main() automatically
#   when included.

using DelimitedFiles
using Printf
using Statistics
using LinearAlgebra
using GLMakie

const BASEDIR = get(ENV, "YZGO_1D_SCAN_DIR",
    raw"C:\Users\vdp\ORNL Dropbox\Daniel Pajerowski\YZGO\CNCS_data\SYM_1d_scans")

outdir_for(base_dir::AbstractString) = joinpath(base_dir, "plots")

struct HistMeta
    label::String
    de_range::Tuple{Float64,Float64}
    plot_ylim_raw::Tuple{Float64,Float64}
end

const Q_ORDER = ["0_1_0", "0p33_0p33_0", "0p5_0_0"]
const FIELD_ORDER_007K = [0.0, 9.0, 14.0]
const CONFIG_ORDER_3P32 = [(0.07, 0.0), (0.07, 9.0), (0.07, 14.0), (20.0, 0.0)]

const QMETA = Dict(
    "0_1_0" => HistMeta("Q = (0, 1, 0)", (0.0, 3.3), (0.0, 0.004)),
    "0p33_0p33_0" => HistMeta("Q = (0.33, 0.33, 0)", (0.0, 4.0), (0.0, 0.004)),
    "0p5_0_0" => HistMeta("Q = (0.5, 0, 0)", (0.0, 4.0), (0.0, 0.004)),
)

struct Scan1D
    qtag::String
    field_T::Float64
    temperature_K::Float64
    Ei_meV::Float64
    filename::String
    intensity::Vector{Float64}
    error::Vector{Float64}
    energy::Vector{Float64}
    col4::Vector{Float64}
    col5::Vector{Float64}
    col6::Vector{Float64}
end

struct BackgroundModel
    qtag::String
    Ei_meV::Float64
    method::String
    energy_grid::Vector{Float64}
    background::Vector{Float64}
    anchor_energy::Vector{Float64}
    anchor_intensity::Vector{Float64}
end

struct PchipInterpolator
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

parsefloat_token(s::AbstractString) = parse(Float64, replace(String(s), "p" => "."))

function parse_filename(path::AbstractString)
    fname = basename(path)
    m = match(r"^yzgo_(\d+p\d+)meV_((?:\d+p\d+)|(?:\d+))K_(\d+(?:p\d+)?)T_Escan_(.+)_SYM\.dat$", fname)
    m === nothing && error("Filename does not match expected YZGO Escan pattern: $fname")

    Ei_meV = parsefloat_token(m.captures[1])
    temperature_K = parsefloat_token(m.captures[2])
    field_T = parsefloat_token(m.captures[3])
    qtag = m.captures[4]
    haskey(QMETA, qtag) || error("No histogram metadata defined for qtag = $qtag")

    return Ei_meV, temperature_K, field_T, qtag
end

function load_scan(path::AbstractString)
    Ei_meV, temperature_K, field_T, qtag = parse_filename(path)

    # File columns are expected to be:
    #   Intensity, Error, DeltaE, coord4, coord5, coord6
    # The coord column order differs between the uploaded 4.65 and 3.32 meV files,
    # so the plotting code deliberately uses the filename qtag for labeling.
    A = readdlm(path, Float64; comments=true, comment_char='#')
    size(A, 2) >= 6 || error("Expected at least 6 columns in $(basename(path)); got $(size(A, 2))")

    return Scan1D(
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

function resolve_base_dir(base_dir::AbstractString)
    if isdir(base_dir)
        return String(base_dir)
    end

    # Handy fallback if this script is placed directly inside the data directory.
    local_dir = @__DIR__
    if isdir(local_dir) && any(f -> endswith(f, ".dat") && occursin("_Escan_", f), readdir(local_dir))
        @warn "BASEDIR was not found; using the script directory instead" base_dir local_dir
        return local_dir
    end

    error("Data directory not found: $base_dir\nEdit BASEDIR, set ENV[\"YZGO_1D_SCAN_DIR\"], or place this script in the data directory.")
end

function load_scans(base_dir::AbstractString=BASEDIR)
    base_dir = resolve_base_dir(base_dir)
    files = filter(readdir(base_dir; join=true)) do f
        fname = basename(f)
        endswith(fname, ".dat") && startswith(fname, "yzgo_") && occursin("_Escan_", fname) && occursin("_SYM.dat", fname)
    end

    scans = Scan1D[]
    for f in sort(files)
        try
            push!(scans, load_scan(f))
        catch err
            @warn "Skipping file that did not parse/load cleanly" file=basename(f) exception=(err, catch_backtrace())
        end
    end

    isempty(scans) && error("No YZGO Escan files were loaded from $base_dir")
    return scans
end

function print_scan_summary(scans::Vector{Scan1D})
    println("\nLoaded YZGO 1D scans")
    println("---------------------")
    for Ei in sort(unique(s.Ei_meV for s in scans))
        println(@sprintf("\nEi = %.2f meV", Ei))
        for qtag in Q_ORDER
            sq = [s for s in scans if isapprox(s.Ei_meV, Ei; atol=1e-8) && s.qtag == qtag]
            isempty(sq) && continue
            println("  ", QMETA[qtag].label)
            for s in sort(sq; by=s -> (s.temperature_K, s.field_T))
                println(@sprintf("    T = %5.2f K, B = %4.1f T: %3d points, E = %.3f to %.3f meV, I = %.4g to %.4g  [%s]",
                    s.temperature_K, s.field_T, length(s.energy), minimum(s.energy), maximum(s.energy),
                    minimum(s.intensity), maximum(s.intensity), s.filename))
            end
        end
    end
end

function scan_exists(scans::Vector{Scan1D}; Ei::Real, T::Real, B::Real, qtag::String)
    return any(s -> isapprox(s.Ei_meV, Ei; atol=1e-8) &&
                    isapprox(s.temperature_K, T; atol=1e-8) &&
                    isapprox(s.field_T, B; atol=1e-8) &&
                    s.qtag == qtag, scans)
end

function get_scan(scans::Vector{Scan1D}; Ei::Real, T::Real, B::Real, qtag::String)
    idx = findfirst(s -> isapprox(s.Ei_meV, Ei; atol=1e-8) &&
                         isapprox(s.temperature_K, T; atol=1e-8) &&
                         isapprox(s.field_T, B; atol=1e-8) &&
                         s.qtag == qtag, scans)
    idx === nothing && error(@sprintf("Missing scan: Ei=%.3f meV, T=%.3f K, B=%.3f T, qtag=%s", Ei, T, B, qtag))
    return scans[idx]
end

function scans_for(scans::Vector{Scan1D}; Ei::Real, qtag::String, T::Union{Nothing,Real}=nothing)
    out = [s for s in scans if isapprox(s.Ei_meV, Ei; atol=1e-8) && s.qtag == qtag]
    if T !== nothing
        out = [s for s in out if isapprox(s.temperature_K, T; atol=1e-8)]
    end
    return sort(out; by=s -> (s.temperature_K, s.field_T))
end

function assert_common_energy_grid(scan_list::Vector{Scan1D}; atol=1e-10)
    isempty(scan_list) && error("No scans supplied for common energy grid check")
    Eref = scan_list[1].energy
    for s in scan_list[2:end]
        length(s.energy) == length(Eref) || error("Energy grid length mismatch for $(s.filename)")
        maximum(abs.(s.energy .- Eref)) <= atol || error("Energy grids are not identical for $(s.filename)")
    end
    return copy(Eref)
end

function sort_xy(xin::AbstractVector{<:Real}, yin::AbstractVector{<:Real})
    length(xin) == length(yin) || error("x and y must have the same length")
    p = sortperm(xin)
    return Float64.(xin[p]), Float64.(yin[p])
end

function gaussian_smooth_xy(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, sigma::Real)
    length(x) == length(y) || error("x and y must have same length")
    sigma <= 0 && return Float64.(y)

    ys = similar(Float64.(y))
    sigma2 = Float64(sigma)^2
    for i in eachindex(x)
        wsum = 0.0
        ysum = 0.0
        xi = Float64(x[i])
        for j in eachindex(x)
            dx = Float64(x[j]) - xi
            w = exp(-0.5 * dx^2 / sigma2)
            wsum += w
            ysum += w * Float64(y[j])
        end
        ys[i] = ysum / wsum
    end
    return ys
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

    h00 = 2.0*t^3 - 3.0*t^2 + 1.0
    h10 = t^3 - 2.0*t^2 + t
    h01 = -2.0*t^3 + 3.0*t^2
    h11 = t^3 - t^2

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
        pchip_interpolator(Es, Is_smooth)(Egrid)
    else
        error("Unknown interpolation_kind=$(interpolation_kind). Use :pchip or :linear.")
    end

    return Float64.(bg), Es, Is_smooth
end

# -------------------------------------------------------------------------
# Ei = 4.65 meV background model, copied/adapted from the starter logic
# -------------------------------------------------------------------------

const MIN_BG_LOW_WINDOW = (0.0, 0.75)
const MIN_BG_HIGH_THRESHOLD = 2.5
const STRUCTURED_FIT_WINDOW = (1.0, 3.0)

const STRUCTURED_RESIDUAL_WINDOWS = Dict(
    "0p33_0p33_0" => (1.675, 2.375),
    "0p5_0_0" => (1.825, 2.425),
)

function ridge_linear_least_squares(X::AbstractMatrix{<:Real},
                                    y::AbstractVector{<:Real},
                                    err::Union{Nothing,AbstractVector{<:Real}}=nothing;
                                    ridge_lambda::Real=0.0,
                                    unpenalized_columns::Vector{Int}=Int[])
    Xf = Float64.(X)
    yf = Float64.(y)
    if err !== nothing
        sigma = max.(Float64.(err), eps(Float64))
        sw = 1.0 ./ sigma
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
        push!(cols, exp.(-Ef ./ Float64(power)))
    elseif model == :inverse_series
        push!(cols, (Ef .+ Float64(energy_offset)) .^ (-1.0))
        push!(cols, (Ef .+ Float64(energy_offset)) .^ (-2.0))
    elseif model == :line
        push!(cols, Ef .- mean(Ef))
    else
        error("Unknown continuum model $(model).")
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
        actual_param = isnan(param) ? 3.0 : Float64(param)
        X = continuum_design_matrix(Ef;
            model=model,
            power=actual_param,
            energy_offset=energy_offset,
            include_linear_tilt=include_linear_tilt,
        )
        coeff = ridge_linear_least_squares(X, If, errfit;
            ridge_lambda=ridge_lambda,
            unpenalized_columns=[1],
        )

        if positive_tail && (model == :power_tail || model == :exp_tail) && coeff[2] < 0
            continue
        end

        pred = X * coeff
        resid = If .- pred
        score = if errfit !== nothing
            sigma = max.(Float64.(errfit), eps(Float64))
            mean((resid ./ sigma).^2)
        else
            mean(resid.^2)
        end

        if score < best_score
            best_score = score
            best_coeff = Float64.(coeff)
            best_param = actual_param
        end
    end

    if isempty(best_coeff)
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

function min_over_scans_background(scan_list::Vector{Scan1D})
    Egrid = assert_common_energy_grid(scan_list)
    Imin = similar(Egrid)
    for i in eachindex(Egrid)
        Imin[i] = minimum(s.intensity[i] for s in scan_list)
    end
    return Egrid, Imin
end

function structured_residual_points(scan0T::Scan1D,
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
    E = scan0T.energy
    I0 = scan0T.intensity
    fit_lo, fit_hi = fit_window
    res_lo, res_hi = residual_window

    peak_mask = (E .>= res_lo) .& (E .<= res_hi)
    fit_mask = (E .>= fit_lo) .& (E .<= fit_hi) .& .!peak_mask

    any(fit_mask) || error("No points found for 0 T continuum fit outside residual window")
    any(peak_mask) || error("No points found inside residual window $(residual_window)")

    Efit = E[fit_mask]
    Ifit = I0[fit_mask]
    errfit = scan0T.error[fit_mask]

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

    @printf("%s 4.65 meV 0T continuum baseline: model=%s, param=%.4g, score=%.4g, coeff=%s\n",
            scan0T.qtag, String(model), param, score, repr(coeff))

    return Eres, residual
end

function make_4p65_background_model(scans::Vector{Scan1D}, qtag::String;
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
    scan_list = [get_scan(scans; Ei=4.65, T=0.07, B=B, qtag=qtag) for B in FIELD_ORDER_007K]
    Egrid, Imin = min_over_scans_background(scan_list)

    lo, hi = MIN_BG_LOW_WINDOW
    anchor_mask = ((Egrid .>= lo) .& (Egrid .<= hi)) .| (Egrid .> MIN_BG_HIGH_THRESHOLD)
    Eraw = collect(Egrid[anchor_mask])
    Iraw = collect(Imin[anchor_mask])

    if haskey(STRUCTURED_RESIDUAL_WINDOWS, qtag)
        scan0T = get_scan(scans; Ei=4.65, T=0.07, B=0.0, qtag=qtag)
        Eres, residual = structured_residual_points(scan0T, STRUCTURED_RESIDUAL_WINDOWS[qtag];
            model=zeroT_baseline_model,
            power_grid=power_grid,
            exp_tau_grid=exp_tau_grid,
            energy_offset=energy_offset,
            ridge_lambda=ridge_lambda,
            include_linear_tilt=include_linear_tilt,
            positive_tail=positive_tail,
            clip_negative_residuals=clip_negative_residuals,
        )

        bg_without_residual, _, _ = make_interpolated_background(Egrid, Eraw, Iraw;
            smooth_sigma_meV=final_smooth_sigma_meV,
            interpolation_kind=final_interp_kind,
        )
        residual_base = linear_interpolate(Egrid, bg_without_residual, Eres)
        residual_abs_background = residual_base .+ residual

        append!(Eraw, Eres)
        append!(Iraw, residual_abs_background)
    end

    bg, Esorted, Ismoothed = make_interpolated_background(Egrid, Eraw, Iraw;
        smooth_sigma_meV=final_smooth_sigma_meV,
        interpolation_kind=final_interp_kind,
    )

    return BackgroundModel(qtag, 4.65, "4.65 meV tail/min-field background", Egrid, bg, Esorted, Ismoothed)
end

# -------------------------------------------------------------------------
# Ei = 3.32 meV background model: point-by-point minimum over configurations
# -------------------------------------------------------------------------

function make_3p32_min_background_model(scans::Vector{Scan1D}, qtag::String;
                                        configs=CONFIG_ORDER_3P32)
    scan_list = [get_scan(scans; Ei=3.32, T=T, B=B, qtag=qtag) for (T, B) in configs]
    Egrid, Imin = min_over_scans_background(scan_list)
    return BackgroundModel(qtag, 3.32, "3.32 meV pointwise minimum over 4 configurations", Egrid, Imin, Egrid, Imin)
end

function subtract_background(scan::Scan1D, bg::BackgroundModel)
    bgvec = if length(scan.energy) == length(bg.energy_grid) && maximum(abs.(scan.energy .- bg.energy_grid)) <= 1e-10
        bg.background
    else
        linear_interpolate(bg.energy_grid, bg.background, scan.energy)
    end

    return Scan1D(
        scan.qtag,
        scan.field_T,
        scan.temperature_K,
        scan.Ei_meV,
        scan.filename,
        scan.intensity .- bgvec,
        copy(scan.error),
        copy(scan.energy),
        copy(scan.col4),
        copy(scan.col5),
        copy(scan.col6),
    )
end

function make_background_models(scans::Vector{Scan1D}; kwargs...)
    models = Dict{Tuple{Float64,String},BackgroundModel}()

    for qtag in Q_ORDER
        if all(B -> scan_exists(scans; Ei=4.65, T=0.07, B=B, qtag=qtag), FIELD_ORDER_007K)
            models[(4.65, qtag)] = make_4p65_background_model(scans, qtag; kwargs...)
        end

        # Only build the 3.32 model if this qtag exists at Ei = 3.32.
        if any(s -> isapprox(s.Ei_meV, 3.32; atol=1e-8) && s.qtag == qtag, scans)
            models[(3.32, qtag)] = make_3p32_min_background_model(scans, qtag)
        end
    end

    return models
end

function make_background_subtracted_scans(scans::Vector{Scan1D}, models::Dict{Tuple{Float64,String},BackgroundModel})
    corrected = Scan1D[]
    for s in scans
        key = (s.Ei_meV, s.qtag)
        if haskey(models, key)
            push!(corrected, subtract_background(s, models[key]))
        end
    end
    return corrected
end

function scan_label(s::Scan1D)
    if isapprox(s.temperature_K, 0.07; atol=1e-8)
        return @sprintf("%.0f T", s.field_T)
    else
        return @sprintf("%.0f K, %.0f T", s.temperature_K, s.field_T)
    end
end

function scan_label_with_ei(s::Scan1D)
    if isapprox(s.temperature_K, 0.07; atol=1e-8)
        return @sprintf("Ei %.2f, %.0f T", s.Ei_meV, s.field_T)
    else
        return @sprintf("Ei %.2f, %.0f K, %.0f T", s.Ei_meV, s.temperature_K, s.field_T)
    end
end

function plot_scan_set!(ax, scan_list::Vector{Scan1D}; show_errorbars::Bool=true)
    for s in sort(scan_list; by=s -> (s.temperature_K, s.field_T))
        ls = isapprox(s.temperature_K, 20.0; atol=1e-8) ? :dash : :solid
        lines!(ax, s.energy, s.intensity; label=scan_label(s), linewidth=2, linestyle=ls)
        if show_errorbars
            errorbars!(ax, s.energy, s.intensity, s.error; whiskerwidth=4)
        end
    end
end

function plot_overlay_scan_set!(ax, scan_list::Vector{Scan1D}; show_errorbars::Bool=true)
    # Linestyle encodes incident energy/configuration:
    #   4.65 meV, 0.07 K : solid
    #   3.32 meV, 0.07 K : dash
    #   3.32 meV, 20 K   : dot
    ordered = sort(scan_list; by=s -> (s.field_T, s.temperature_K, s.Ei_meV))
    for s in ordered
        ls = if isapprox(s.Ei_meV, 4.65; atol=1e-8)
            :solid
        elseif isapprox(s.temperature_K, 20.0; atol=1e-8)
            :dot
        else
            :dash
        end
        lw = isapprox(s.Ei_meV, 4.65; atol=1e-8) ? 2.4 : 2.0
        lines!(ax, s.energy, s.intensity; label=scan_label_with_ei(s), linewidth=lw, linestyle=ls)
        if show_errorbars
            errorbars!(ax, s.energy, s.intensity, s.error; whiskerwidth=4)
        end
    end
end

function plot_bgsub_ei_comparison(bgsub_scans::Vector{Scan1D};
                                  outdir::AbstractString=outdir_for(resolve_base_dir(BASEDIR)),
                                  save_png::Bool=true,
                                  filename::String="YZGO_1d_Escan_4p65_3p32_bgsub_comparison.png",
                                  bgsub_ylim::Tuple{Float64,Float64}=(-0.0005, 0.0030),
                                  show_errorbars::Bool=true)
    mkpath(outdir)
    fig = Figure(size=(1500, 1180))

    for (icol, qtag) in enumerate(Q_ORDER)
        meta = QMETA[qtag]

        # Rows 1 and 2: separate Ei plots.
        for (irow, Ei) in enumerate((4.65, 3.32))
            ax = Axis(fig[irow, icol];
                title = irow == 1 ? meta.label : "",
                xlabel = "",
                ylabel = icol == 1 ? @sprintf("Ei = %.2f meV\nI - background", Ei) : "",
                limits = (meta.de_range[1], meta.de_range[2], bgsub_ylim[1], bgsub_ylim[2]),
            )

            scan_list = scans_for(bgsub_scans; Ei=Ei, qtag=qtag)
            if isempty(scan_list)
                text!(ax, 0.5, 0.5; text="not measured", space=:relative, align=(:center, :center), fontsize=22)
                hidedecorations!(ax; grid=false)
                hidespines!(ax)
                continue
            end

            plot_scan_set!(ax, scan_list; show_errorbars=show_errorbars)
            hlines!(ax, [0.0]; linestyle=:dash, linewidth=1)
            axislegend(ax; position=:rt, framevisible=false)
        end

        # Row 3: direct overlay of all available background-subtracted scans for this Q.
        ax_overlay = Axis(fig[3, icol];
            xlabel = "DeltaE (meV)",
            ylabel = icol == 1 ? "overlay\nI - background" : "",
            limits = (meta.de_range[1], meta.de_range[2], bgsub_ylim[1], bgsub_ylim[2]),
        )

        overlay_list = vcat(
            scans_for(bgsub_scans; Ei=4.65, qtag=qtag),
            scans_for(bgsub_scans; Ei=3.32, qtag=qtag),
        )

        if isempty(overlay_list)
            text!(ax_overlay, 0.5, 0.5; text="not measured", space=:relative, align=(:center, :center), fontsize=22)
            hidedecorations!(ax_overlay; grid=false)
            hidespines!(ax_overlay)
            continue
        end

        plot_overlay_scan_set!(ax_overlay, overlay_list; show_errorbars=show_errorbars)
        hlines!(ax_overlay, [0.0]; linestyle=:dash, linewidth=1)

        if isempty(scans_for(bgsub_scans; Ei=3.32, qtag=qtag)) && !isempty(scans_for(bgsub_scans; Ei=4.65, qtag=qtag))
            text!(ax_overlay, 0.5, 0.88;
                text="Ei = 3.32 meV not measured here",
                space=:relative,
                align=(:center, :center),
                fontsize=16,
            )
        end

        axislegend(ax_overlay; position=:rt, framevisible=false)
    end

    Label(fig[0, :], "YZGO 1D scans: background-subtracted Ei = 4.65 meV, Ei = 3.32 meV, and direct overlay", fontsize=22)

    if save_png
        save(joinpath(outdir, filename), fig)
        println("Saved ", joinpath(outdir, filename))
    end

    return fig
end

function plot_background_models(bg_models::Dict{Tuple{Float64,String},BackgroundModel};
                                outdir::AbstractString=outdir_for(resolve_base_dir(BASEDIR)),
                                save_png::Bool=true,
                                filename::String="YZGO_1d_Escan_background_models_4p65_3p32.png")
    mkpath(outdir)
    fig = Figure(size=(1500, 860))

    for (icol, qtag) in enumerate(Q_ORDER)
        meta = QMETA[qtag]
        for (irow, Ei) in enumerate((4.65, 3.32))
            ax = Axis(fig[irow, icol];
                title = irow == 1 ? meta.label : "",
                xlabel = irow == 2 ? "DeltaE (meV)" : "",
                ylabel = icol == 1 ? @sprintf("Ei = %.2f meV\nbackground", Ei) : "",
            )
            xlims!(ax, meta.de_range[1], meta.de_range[2])

            key = (Ei, qtag)
            if !haskey(bg_models, key)
                text!(ax, 0.5, 0.5; text="not measured", space=:relative, align=(:center, :center), fontsize=22)
                hidedecorations!(ax; grid=false)
                hidespines!(ax)
                continue
            end

            bg = bg_models[key]
            scatter!(ax, bg.anchor_energy, bg.anchor_intensity; label="anchors", markersize=6)
            lines!(ax, bg.energy_grid, bg.background; label=bg.method, linewidth=3)
            axislegend(ax; position=:rt, framevisible=false)
        end
    end

    Label(fig[0, :], "YZGO 1D background models", fontsize=22)

    if save_png
        save(joinpath(outdir, filename), fig)
        println("Saved ", joinpath(outdir, filename))
    end

    return fig
end

function main(; base_dir::AbstractString=BASEDIR,
                outdir::AbstractString=outdir_for(resolve_base_dir(base_dir)),
                save_png::Bool=true,
                show_errorbars::Bool=true,
                bgsub_ylim::Tuple{Float64,Float64}=(-0.0005, 0.0030),
                make_diagnostic_plot::Bool=true,
                display_figures::Bool=false,
                background_kwargs...)
    base_dir = resolve_base_dir(base_dir)
    mkpath(outdir)

    scans = load_scans(base_dir)
    print_scan_summary(scans)

    bg_models = make_background_models(scans; background_kwargs...)
    bgsub_scans = make_background_subtracted_scans(scans, bg_models)

    fig = plot_bgsub_ei_comparison(bgsub_scans;
        outdir=outdir,
        save_png=save_png,
        bgsub_ylim=bgsub_ylim,
        show_errorbars=show_errorbars,
    )

    if make_diagnostic_plot
        plot_background_models(bg_models; outdir=outdir, save_png=save_png)
    end

    if display_figures
        display(fig)
    end

    return scans, bgsub_scans, bg_models, fig
end

# Do not run automatically when included by repo-native driver scripts.
#
# Running directly from the command line is still supported:
#
#     julia scripts/legacy/yzgo_plot_1d_scans_4p65_3p32_compare_legacy.jl
#
if abspath(PROGRAM_FILE) == @__FILE__
    GLMakie.activate!()
    main(display_figures=true)
end
