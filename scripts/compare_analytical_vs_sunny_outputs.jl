# scripts/compare_analytical_vs_sunny_outputs.jl
#
# Post-processing comparison of already-saved analytical and Sunny.jl model
# outputs.  This script does not recompute either model; it only reads CSV files
# and makes overplot / map comparison figures.
#
# Run from the repo root with:
#
#   julia --project=. scripts/compare_analytical_vs_sunny_outputs.jl
#
# Optional first argument:
#
#   julia --project=. scripts/compare_analytical_vs_sunny_outputs.jl results/fits/cofit_9T14T_shared_fraction

using Printf
using Statistics
using CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const DEFAULT_ANALYTICAL_1D_DIR = joinpath(REPO_ROOT, "results", "fits", "cofit_9T14T_shared_fraction")
const ANALYTICAL_1D_DIR = length(ARGS) >= 1 ? normpath(isabspath(ARGS[1]) ? ARGS[1] : joinpath(REPO_ROOT, ARGS[1])) : DEFAULT_ANALYTICAL_1D_DIR

const SUNNY_TABLE_DIR = joinpath(REPO_ROOT, "results", "feature_tables", "sunny_validation")
const ANALYTICAL_2D_TABLE_DIR = joinpath(REPO_ROOT, "results", "feature_tables", "analytical_2d_model")
const OUTFIG_DIR = joinpath(REPO_ROOT, "results", "figures", "model_comparison")
const OUTTABLE_DIR = joinpath(REPO_ROOT, "results", "feature_tables", "model_comparison")

const FIELDS_T = [9.0, 14.0]
const QTAGS = ["0_1_0", "0p33_0p33_0", "0p5_0_0"]
const QLABELS = Dict(
    "0_1_0" => "Γ / (0,1,0)",
    "0p33_0p33_0" => "K / (1/3,1/3,0)",
    "0p5_0_0" => "M / (1/2,0,0)",
)

# -----------------------------------------------------------------------------
# Minimal CSV utilities
# -----------------------------------------------------------------------------

function csv_cell(x)
    if x === missing || x === nothing
        return ""
    elseif x isa AbstractFloat
        if isnan(x)
            return "NaN"
        elseif isinf(x)
            return x > 0 ? "Inf" : "-Inf"
        else
            return @sprintf("%.12g", x)
        end
    elseif x isa Integer
        return string(x)
    elseif x isa Bool
        return string(x)
    else
        s = string(x)
        if occursin(',', s) || occursin('"', s) || occursin('\n', s)
            return "\"" * replace(s, "\"" => "\"\"") * "\""
        else
            return s
        end
    end
end

function write_csv(path::AbstractString, header::Vector{String}, cols::Vector)
    mkpath(dirname(path))
    isempty(cols) && error("No columns supplied for $path")
    n = length(cols[1])
    for (j, c) in enumerate(cols)
        length(c) == n || error("Column $(header[j]) has length $(length(c)); expected $n")
    end
    open(path, "w") do io
        println(io, join(header, ","))
        for i in 1:n
            println(io, join((csv_cell(c[i]) for c in cols), ","))
        end
    end
    return path
end

function read_csv_table(path::AbstractString)
    isfile(path) || error("Could not find CSV file: $path")
    lines = readlines(path)
    isempty(lines) && error("CSV file is empty: $path")
    header = strip.(split(strip(lines[1]), ","))
    cols = Dict(h => String[] for h in header)
    for raw in lines[2:end]
        isempty(strip(raw)) && continue
        parts = strip.(split(raw, ","; keepempty=true))
        if length(parts) < length(header)
            append!(parts, fill("", length(header) - length(parts)))
        elseif length(parts) > length(header)
            parts = parts[1:length(header)]
        end
        for (h, v) in zip(header, parts)
            push!(cols[h], v)
        end
    end
    return (; path, header, cols)
end

nrows(tbl) = isempty(tbl.header) ? 0 : length(tbl.cols[tbl.header[1]])
hascol(tbl, name::AbstractString) = haskey(tbl.cols, name)

function parse_float_or_nan(s)
    try
        return parse(Float64, strip(String(s)))
    catch
        return NaN
    end
end

function colnum(tbl, name::AbstractString; default::Union{Nothing,Real}=nothing)
    if hascol(tbl, name)
        return [parse_float_or_nan(s) for s in tbl.cols[name]]
    elseif default !== nothing
        return fill(Float64(default), nrows(tbl))
    else
        error("Missing column '$name' in $(tbl.path). Available columns: $(join(tbl.header, ", "))")
    end
end

function colstr(tbl, name::AbstractString; default::Union{Nothing,String}=nothing)
    if hascol(tbl, name)
        return tbl.cols[name]
    elseif default !== nothing
        return fill(default, nrows(tbl))
    else
        error("Missing column '$name' in $(tbl.path). Available columns: $(join(tbl.header, ", "))")
    end
end

function mask_field_qtag(tbl; field_T::Real, qtag::AbstractString)
    n = nrows(tbl)
    m = trues(n)
    if hascol(tbl, "field_T")
        f = colnum(tbl, "field_T")
        m .&= abs.(f .- Float64(field_T)) .< 1e-6
    end
    if hascol(tbl, "qtag")
        q = colstr(tbl, "qtag")
        m .&= (q .== String(qtag))
    end
    return m
end

function safe_sort_xy(x, ys...)
    p = sortperm(x)
    return (x[p], (y[p] for y in ys)...)
end

function first_existing(paths::Vector{String})
    for p in paths
        isfile(p) && return p
    end
    return nothing
end

field_tag(B::Real) = @sprintf("%gT", Float64(B))

# -----------------------------------------------------------------------------
# 1D comparison
# -----------------------------------------------------------------------------

function analytical_1d_paths(dir::AbstractString)
    model_path = first_existing([
        joinpath(dir, "YZGO_neutron_magnetization_shared_fraction_cofit_neutron_scaled_model.csv"),
        joinpath(dir, "YZGO_1d_two_kernel_scaled_model.csv"),
    ])
    fit_path = first_existing([
        joinpath(dir, "YZGO_neutron_magnetization_shared_fraction_cofit_neutron_fit_points.csv"),
        joinpath(dir, "YZGO_1d_two_kernel_fit_points.csv"),
    ])
    return (; model_path, fit_path)
end

function sunny_1d_path(qtag::AbstractString, B::Real)
    return joinpath(SUNNY_TABLE_DIR, @sprintf("sunny_kpm_1d_%s_%gT_vs_exp.csv", qtag, Float64(B)))
end

function write_1d_long_csv(analytical_model, analytical_fit, sunny_tables)
    model_name = String[]
    field_T = Float64[]
    qtag_col = String[]
    energy_meV = Float64[]
    I_exp = Float64[]
    Ierr_exp = Float64[]
    I_total = Float64[]
    I_disp = Float64[]
    I_flat = Float64[]
    source_csv = String[]

    if analytical_model !== nothing
        for B in FIELDS_T, qtag in QTAGS
            m = mask_field_qtag(analytical_model; field_T=B, qtag=qtag)
            if any(m)
                E = colnum(analytical_model, "energy_meV")[m]
                It = colnum(analytical_model, "model_total_scaled")[m]
                Id = colnum(analytical_model, "model_dispersive_scaled")[m]
                If = colnum(analytical_model, "model_nondispersive_scaled")[m]
                for i in eachindex(E)
                    push!(model_name, "analytical")
                    push!(field_T, B); push!(qtag_col, qtag); push!(energy_meV, E[i])
                    push!(I_exp, NaN); push!(Ierr_exp, NaN)
                    push!(I_total, It[i]); push!(I_disp, Id[i]); push!(I_flat, If[i])
                    push!(source_csv, analytical_model.path)
                end
            end
        end
    end

    for ((B, qtag), tbl) in sunny_tables
        E = colnum(tbl, "energy_meV")
        It = colnum(tbl, "I_total_scaled")
        Id = colnum(tbl, "I_disp_scaled")
        If = colnum(tbl, "I_flat_scaled")
        Ie = hascol(tbl, "I_exp") ? colnum(tbl, "I_exp") : fill(NaN, nrows(tbl))
        Ier = hascol(tbl, "Ierr_exp") ? colnum(tbl, "Ierr_exp") : fill(NaN, nrows(tbl))
        for i in eachindex(E)
            push!(model_name, "sunny")
            push!(field_T, B); push!(qtag_col, qtag); push!(energy_meV, E[i])
            push!(I_exp, Ie[i]); push!(Ierr_exp, Ier[i])
            push!(I_total, It[i]); push!(I_disp, Id[i]); push!(I_flat, If[i])
            push!(source_csv, tbl.path)
        end
    end

    path = joinpath(OUTTABLE_DIR, "analytical_vs_sunny_1d_long.csv")
    write_csv(path,
        ["model", "field_T", "qtag", "energy_meV", "I_exp", "Ierr_exp", "I_total", "I_disp", "I_flat", "source_csv"],
        [model_name, field_T, qtag_col, energy_meV, I_exp, Ierr_exp, I_total, I_disp, I_flat, source_csv],
    )
    return path
end

function plot_1d_comparison()
    mkpath(OUTFIG_DIR); mkpath(OUTTABLE_DIR)
    paths = analytical_1d_paths(ANALYTICAL_1D_DIR)
    analytical_model = paths.model_path === nothing ? nothing : read_csv_table(paths.model_path)
    analytical_fit = paths.fit_path === nothing ? nothing : read_csv_table(paths.fit_path)

    sunny_tables = Dict{Tuple{Float64,String},Any}()
    for B in FIELDS_T, qtag in QTAGS
        p = sunny_1d_path(qtag, B)
        if isfile(p)
            sunny_tables[(B, qtag)] = read_csv_table(p)
        else
            @warn "Missing Sunny 1D CSV" p
        end
    end

    if analytical_model === nothing
        @warn "No analytical 1D model CSV found" dir=ANALYTICAL_1D_DIR
    end
    if isempty(sunny_tables)
        @warn "No Sunny 1D CSVs found" dir=SUNNY_TABLE_DIR
    end

    long_path = write_1d_long_csv(analytical_model, analytical_fit, sunny_tables)

    fig = Figure(size=(1500, 760), fontsize=14)
    Label(fig[0, :], "Analytical vs Sunny KPM 1D model outputs", fontsize=22, font=:bold)

    for (irow, B) in enumerate(FIELDS_T), (icol, qtag) in enumerate(QTAGS)
        ax = Axis(fig[irow, icol],
            title=@sprintf("%g T, %s", B, get(QLABELS, qtag, qtag)),
            xlabel="Energy transfer ΔE (meV)",
            ylabel=icol == 1 ? "Intensity (arb.)" : "",
        )

        # Experiment from Sunny-vs-exp CSV if available. This is the full scan grid.
        stbl = get(sunny_tables, (B, qtag), nothing)
        if stbl !== nothing && hascol(stbl, "I_exp")
            E = colnum(stbl, "energy_meV")
            y = colnum(stbl, "I_exp")
            yerr = hascol(stbl, "Ierr_exp") ? colnum(stbl, "Ierr_exp") : fill(NaN, length(E))
            scatter!(ax, E, y; markersize=5, label="experiment")
            if any(isfinite, yerr)
                errorbars!(ax, E, y, yerr; whiskerwidth=3)
            end
        elseif analytical_fit !== nothing
            m = mask_field_qtag(analytical_fit; field_T=B, qtag=qtag)
            if any(m)
                E = colnum(analytical_fit, "energy_meV")[m]
                y = colnum(analytical_fit, "data_intensity")[m]
                yerr = hascol(analytical_fit, "data_error") ? colnum(analytical_fit, "data_error")[m] : fill(NaN, length(E))
                scatter!(ax, E, y; markersize=5, label="experiment fit pts")
                any(isfinite, yerr) && errorbars!(ax, E, y, yerr; whiskerwidth=3)
            end
        end

        # Analytical model.
        if analytical_model !== nothing
            m = mask_field_qtag(analytical_model; field_T=B, qtag=qtag)
            if any(m)
                E = colnum(analytical_model, "energy_meV")[m]
                It = colnum(analytical_model, "model_total_scaled")[m]
                Id = colnum(analytical_model, "model_dispersive_scaled")[m]
                If = colnum(analytical_model, "model_nondispersive_scaled")[m]
                E, It, Id, If = safe_sort_xy(E, It, Id, If)
                lines!(ax, E, It; linewidth=3, label="analytical total")
                lines!(ax, E, Id; linewidth=2, linestyle=:dash, label="analytical disp")
                lines!(ax, E, If; linewidth=2, linestyle=:dot, label="analytical flat")
            end
        end

        # Sunny model.
        if stbl !== nothing
            E = colnum(stbl, "energy_meV")
            It = colnum(stbl, "I_total_scaled")
            Id = hascol(stbl, "I_disp_scaled") ? colnum(stbl, "I_disp_scaled") : fill(NaN, length(E))
            If = hascol(stbl, "I_flat_scaled") ? colnum(stbl, "I_flat_scaled") : fill(NaN, length(E))
            E, It, Id, If = safe_sort_xy(E, It, Id, If)
            lines!(ax, E, It; linewidth=3, linestyle=:dashdot, label="Sunny total")
            lines!(ax, E, Id; linewidth=2, linestyle=:dash, label="Sunny disp")
            lines!(ax, E, If; linewidth=2, linestyle=:dot, label="Sunny flat")
        end

        xlims!(ax, 0.0, 3.3)
        ylims!(ax, -5e-4, 3e-3)
        if irow == 1 && icol == length(QTAGS)
            axislegend(ax; position=:rt, framevisible=false, labelsize=9)
        end
    end

    fig_path = joinpath(OUTFIG_DIR, "analytical_vs_sunny_1d_9T14T.png")
    save(fig_path, fig)
    return (; fig_path, long_path, analytical_model_path=paths.model_path, analytical_fit_path=paths.fit_path, sunny_count=length(sunny_tables))
end

# -----------------------------------------------------------------------------
# 2D comparison
# -----------------------------------------------------------------------------

function pivot_long_table(tbl; xcol="path_coordinate", ecol="energy_meV", zcol="I_total_scaled")
    x = colnum(tbl, xcol)
    e = colnum(tbl, ecol)
    z = colnum(tbl, zcol)
    xs = sort(unique(x[isfinite.(x)]))
    es = sort(unique(e[isfinite.(e)]))
    Z = fill(NaN, length(xs), length(es))
    xidx = Dict(v => i for (i, v) in enumerate(xs))
    eidx = Dict(v => i for (i, v) in enumerate(es))
    for i in eachindex(z)
        if isfinite(x[i]) && isfinite(e[i]) && haskey(xidx, x[i]) && haskey(eidx, e[i])
            Z[xidx[x[i]], eidx[e[i]]] = z[i]
        end
    end
    return (; x=xs, e=es, z=Z)
end

function analytical_2d_path(B::Real; leg::Integer=1)
    return joinpath(ANALYTICAL_2D_TABLE_DIR, @sprintf("analytical_2d_%gT_leg%d.csv", Float64(B), leg))
end

function sunny_2d_path(B::Real; leg::Integer=1)
    return joinpath(SUNNY_TABLE_DIR, @sprintf("sunny_kpm_2d_data_model_%gT_leg%d.csv", Float64(B), leg))
end

function robust_colorrange(arrays; qhi=0.995)
    vals = Float64[]
    for A in arrays
        append!(vals, vec(A)[isfinite.(vec(A))])
    end
    isempty(vals) && return (0.0, 1.0)
    lo = min(0.0, quantile(vals, 0.01))
    hi = quantile(vals, qhi)
    if !(isfinite(hi) && hi > lo)
        hi = maximum(vals)
    end
    if !(isfinite(hi) && hi > lo)
        hi = lo + 1.0
    end
    return (lo, hi)
end

function plot_2d_comparison(; leg::Integer=1)
    mkpath(OUTFIG_DIR); mkpath(OUTTABLE_DIR)
    analytical = Dict{Float64,Any}()
    sunny = Dict{Float64,Any}()
    for B in FIELDS_T
        ap = analytical_2d_path(B; leg)
        sp = sunny_2d_path(B; leg)
        if isfile(ap)
            analytical[B] = pivot_long_table(read_csv_table(ap); zcol="I_model_scaled")
        else
            @warn "Missing analytical 2D CSV. Run scripts/export_analytical_2d_model_csv.jl first." ap
        end
        if isfile(sp)
            sunny[B] = pivot_long_table(read_csv_table(sp); zcol="I_total_scaled")
        else
            @warn "Missing Sunny 2D CSV" sp
        end
    end

    if isempty(analytical) || isempty(sunny)
        return (; fig_path=nothing, analytical_count=length(analytical), sunny_count=length(sunny))
    end

    fields_have = [B for B in FIELDS_T if haskey(analytical, B) && haskey(sunny, B)]
    isempty(fields_have) && return (; fig_path=nothing, analytical_count=length(analytical), sunny_count=length(sunny))

    model_cr = robust_colorrange(vcat([analytical[B].z for B in fields_have], [sunny[B].z for B in fields_have]); qhi=0.995)

    fig = Figure(size=(1200, 900), fontsize=15)
    Label(fig[0, 1:(length(fields_have)+1)], @sprintf("Analytical vs Sunny KPM 2D models, leg %d", leg), fontsize=22, font=:bold, tellwidth=false)
    hm = nothing
    for (icol, B) in enumerate(fields_have)
        a = analytical[B]
        s = sunny[B]
        axa = Axis(fig[1, icol], title=@sprintf("%g T", B), ylabel=icol == 1 ? "Analytical\nΔE (meV)" : "", xlabel="")
        hm = heatmap!(axa, a.x, a.e, a.z; colormap=:viridis, colorrange=model_cr, nan_color=:lightgray)
        ylims!(axa, 0.2, 3.2)
        vlines!(axa, [-1/3, 0.0, 1/3, 2/3, 1.0]; color=(:white, 0.45), linewidth=1)

        axs = Axis(fig[2, icol], title="", ylabel=icol == 1 ? "Sunny KPM\nΔE (meV)" : "", xlabel="Path coordinate (rlu)")
        hm = heatmap!(axs, s.x, s.e, s.z; colormap=:viridis, colorrange=model_cr, nan_color=:lightgray)
        ylims!(axs, 0.2, 3.2)
        vlines!(axs, [-1/3, 0.0, 1/3, 2/3, 1.0]; color=(:white, 0.45), linewidth=1)
    end
    hm !== nothing && Colorbar(fig[1:2, length(fields_have)+1], hm; label="Scaled intensity")
    fig_path = joinpath(OUTFIG_DIR, @sprintf("analytical_vs_sunny_2d_leg%d.png", leg))
    save(fig_path, fig)
    return (; fig_path, analytical_count=length(analytical), sunny_count=length(sunny))
end

function main()
    println("Analytical vs Sunny output comparison")
    println("-------------------------------------")
    println("Analytical 1D directory: ", ANALYTICAL_1D_DIR)
    println("Analytical 2D directory: ", ANALYTICAL_2D_TABLE_DIR)
    println("Sunny directory:         ", SUNNY_TABLE_DIR)
    println("Output figures:          ", OUTFIG_DIR)
    println("Output tables:           ", OUTTABLE_DIR)
    println()

    r1 = plot_1d_comparison()
    println("1D comparison figure: ", r1.fig_path)
    println("1D long CSV:          ", r1.long_path)
    println("Analytical 1D model:  ", r1.analytical_model_path)
    println("Sunny 1D file count:  ", r1.sunny_count)
    println()

    r2 = plot_2d_comparison(; leg=1)
    if r2.fig_path === nothing
        println("2D comparison figure not made. Missing analytical and/or Sunny 2D CSVs.")
        println("Run: julia --project=. scripts/export_analytical_2d_model_csv.jl")
        println("and: julia --project=. scripts/sunny_plot_kpm_2d.jl")
    else
        println("2D comparison figure: ", r2.fig_path)
    end
end

main()
