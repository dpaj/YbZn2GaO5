# scripts/export_analytical_2d_model_csv.jl
#
# Export the analytical fixed-parameter 2D data/model calculation to long CSV
# files so it can be compared directly with Sunny KPM outputs.  This script
# reuses the existing legacy analytical 2D plot implementation and writes one
# CSV per field/leg.
#
# Run from repo root:
#
#   julia --project=. scripts/export_analytical_2d_model_csv.jl

using Printf

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

const PLOT_2D_CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "plot_2d_controls.toml")
const BEST_FIT_PATH = joinpath(REPO_ROOT, "configs", "best_fit_parameters.toml")
const plot_controls = load_toml_config(PLOT_2D_CONTROLS_PATH)

const NEUTRON_2D_DIR = joinpath(REPO_ROOT, splitpath(plot_controls["data"]["neutron_2d_subdir"])...)
const OUTFIG_DIR = joinpath(REPO_ROOT, splitpath(plot_controls["output"]["figure_subdir"])...)
const OUTTABLE_DIR = joinpath(REPO_ROOT, "results", "feature_tables", "analytical_2d_model")

mkpath(OUTFIG_DIR)
mkpath(OUTTABLE_DIR)

# The legacy implementation reads these environment variables at include time.
ENV["YZGO_DATA_DIR"] = NEUTRON_2D_DIR
ENV["YZGO_OUT_DIR"] = OUTFIG_DIR
ENV["MAKIE_BACKEND"] = get(plot_controls["makie"], "backend", "CairoMakie")
ENV["YZGO_BEST_FIT_PARAMETERS_TOML"] = BEST_FIT_PATH

include(joinpath(REPO_ROOT, "scripts", "legacy", "plot_yzgo_2d_data_vs_model_legacy.jl"))

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
    elseif x isa Integer || x isa Bool
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

function maybe_get_matrix_value(A, i, j)
    (i <= size(A, 1) && j <= size(A, 2)) ? A[i, j] : NaN
end

function export_analytical_2d_csvs(result)
    inventory_field_T = Float64[]
    inventory_leg = Int[]
    inventory_path = String[]
    inventory_scale = Float64[]
    inventory_nx = Int[]
    inventory_ne = Int[]

    for key in sort(collect(keys(result.models)))
        leg, field_T = key
        model = result.models[key]
        zscaled = result.scaled_models[key]
        scan = haskey(result.scans, key) ? result.scans[key] : nothing
        scale = result.scale_by_key[key]

        nx = length(model.x)
        ne = length(model.e)
        rows = nx * ne

        field_col = fill(Float64(field_T), rows)
        leg_col = fill(Int(leg), rows)
        path_coordinate = Float64[]
        energy_meV = Float64[]
        I_exp = Float64[]
        I_model_scaled = Float64[]
        I_model_total_unscaled = Float64[]
        I_model_disp_unscaled = Float64[]
        I_model_flat_unweighted = Float64[]
        model_scale = fill(Float64(scale), rows)
        r2 = fill(Float64(model.r2), rows)
        data_grid_matched = fill(scan !== nothing && size(scan.z) == size(model.z), rows)

        for ix in 1:nx, ie in 1:ne
            push!(path_coordinate, model.x[ix])
            push!(energy_meV, model.e[ie])
            push!(I_exp, scan === nothing ? NaN : maybe_get_matrix_value(scan.z, ix, ie))
            push!(I_model_scaled, zscaled[ix, ie])
            push!(I_model_total_unscaled, model.z[ix, ie])
            push!(I_model_disp_unscaled, model.disp.intensity[ix, ie])
            push!(I_model_flat_unweighted, model.flat.intensity[ix, ie])
        end

        ftag = @sprintf("%gT", Float64(field_T))
        path = joinpath(OUTTABLE_DIR, "analytical_2d_$(ftag)_leg$(leg).csv")
        write_csv(path,
            ["field_T", "leg", "path_coordinate", "energy_meV", "I_exp", "I_model_scaled", "I_model_total_unscaled", "I_model_disp_unscaled", "I_model_flat_unweighted", "model_scale", "r2_shared", "data_grid_matched"],
            [field_col, leg_col, path_coordinate, energy_meV, I_exp, I_model_scaled, I_model_total_unscaled, I_model_disp_unscaled, I_model_flat_unweighted, model_scale, r2, data_grid_matched],
        )

        push!(inventory_field_T, Float64(field_T))
        push!(inventory_leg, Int(leg))
        push!(inventory_path, path)
        push!(inventory_scale, Float64(scale))
        push!(inventory_nx, nx)
        push!(inventory_ne, ne)
        println("Wrote: ", path)
    end

    inv_path = joinpath(OUTTABLE_DIR, "analytical_2d_inventory.csv")
    write_csv(inv_path,
        ["field_T", "leg", "csv_path", "model_scale", "n_path", "n_energy"],
        [inventory_field_T, inventory_leg, inventory_path, inventory_scale, inventory_nx, inventory_ne],
    )
    println("Wrote inventory: ", inv_path)
    return inv_path
end

println("Running existing analytical 2D calculation...")
println("  data dir:  ", NEUTRON_2D_DIR)
println("  fig dir:   ", OUTFIG_DIR)
println("  table dir: ", OUTTABLE_DIR)
println()

result = run_yzgo_2d_prelim_model_comparison()
export_analytical_2d_csvs(result)
