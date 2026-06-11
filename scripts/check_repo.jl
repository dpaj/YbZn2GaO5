# scripts/check_repo.jl
#
# Lightweight repository health check for the YbZn2GaO5 analysis repo.
#
# Run from the repo root with:
#
#     julia --project=. scripts/check_repo.jl
#
# This intentionally runs only cheap/no-optimization checks.

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

println()
println("YbZn2GaO5 repo health check")
println("===========================")
println()
println("Repo root:")
println(REPO_ROOT)
println()

# ---------------------------------------------------------------------------
# Check utility module and configs
# ---------------------------------------------------------------------------

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

const BEST_FIT_PATH = joinpath(REPO_ROOT, "configs", "best_fit_parameters.toml")
const COFIT_SMOKE_CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "cofit_controls_smoke.toml")
const PLOT_2D_CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "plot_2d_controls.toml")
const BACKGROUND_COMPARE_CONTROLS_PATH = joinpath(
    REPO_ROOT,
    "configs",
    "background_compare_controls.toml",
)
const SUNNY_VALIDATION_CONTROLS_PATH = joinpath(
    REPO_ROOT,
    "configs",
    "sunny_validation_controls.toml",
)

println("Checking config files...")
for path in (BEST_FIT_PATH, COFIT_SMOKE_CONTROLS_PATH, PLOT_2D_CONTROLS_PATH, BACKGROUND_COMPARE_CONTROLS_PATH, SUNNY_VALIDATION_CONTROLS_PATH)
    if !isfile(path)
        error("Missing required config file: $path")
    end
    println("  found: ", relpath(path, REPO_ROOT))
end
println()

println("Loading canonical model parameters...")
params = load_canonical_model_parameters(BEST_FIT_PATH)
print_canonical_model_parameters(params)

if !hasproperty(params, :magnetization_global_scale)
    error("Canonical parameters did not load magnetization_global_scale")
end
println("  magnetization_global_scale: ", params.magnetization_global_scale)
println()

println("Loading co-fit smoke controls...")
smoke_controls = load_cofit_controls(COFIT_SMOKE_CONTROLS_PATH)
println("  fit name: ", smoke_controls["output"]["fit_name"])
println("  run optimization: ", smoke_controls["optimization"]["run_optimization"])
println()

println("Loading 2D plotting controls...")
plot_controls = load_toml_config(PLOT_2D_CONTROLS_PATH)
println("  2D data subdir: ", plot_controls["data"]["neutron_2d_subdir"])
println("  output subdir: ", plot_controls["output"]["figure_subdir"])
println()

println("Loading Sunny validation controls...")
sunny_controls = load_toml_config(SUNNY_VALIDATION_CONTROLS_PATH)
println("  table subdir: ", sunny_controls["paths"]["table_subdir"])
println("  figure subdir: ", sunny_controls["paths"]["figure_subdir"])
println("  use magnetization scale from best fit: ",
        sunny_controls["magnetization"]["use_magnetization_global_scale_from_best_fit"])
println()

# ---------------------------------------------------------------------------
# Check important data directories/files
# ---------------------------------------------------------------------------

println("Checking data inputs...")

required_paths = [
    joinpath(REPO_ROOT, "data", "neutron", "CNCS_1d_scans"),
    joinpath(REPO_ROOT, "data", "neutron", "CNCS_2d_scans"),
    joinpath(REPO_ROOT, "data", "magnetization", "YZGO_MvB_black_curve_digitized_visible.csv"),
    joinpath(REPO_ROOT, "data", "cif", "ICSD_CollCode138765.cif"),
]

for path in required_paths
    if !ispath(path)
        error("Missing required data path: $path")
    end
    println("  found: ", relpath(path, REPO_ROOT))
end
println()

# ---------------------------------------------------------------------------
# Summarize neutron data counts
# ---------------------------------------------------------------------------

one_d_dir = joinpath(REPO_ROOT, "data", "neutron", "CNCS_1d_scans")
two_d_dir = joinpath(REPO_ROOT, "data", "neutron", "CNCS_2d_scans")

one_d_files = filter(f -> endswith(lowercase(f), ".dat"), readdir(one_d_dir))
two_d_files = filter(f -> endswith(lowercase(f), ".dat"), readdir(two_d_dir))

println("Data file counts:")
println("  1D neutron .dat files: ", length(one_d_files))
println("  2D neutron .dat files: ", length(two_d_files))
println()

println("Repo health check completed successfully.")
println()
println("Next optional checks:")
println("  julia --project=. scripts/run_cofit_9T14T_smoke.jl")
println("  julia --project=. scripts/plot_2d_data_vs_model.jl")
println("  julia --project=. scripts/compare_1d_4p65_3p32_backgrounds.jl")