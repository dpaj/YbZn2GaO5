# scripts/plot_2d_data_vs_model.jl
#
# Repo-native driver for the YbZn2GaO5 2D data-versus-model plot.
#
# This driver reads:
#
#   configs/plot_2d_controls.toml     2D plotting paths/backend controls
#   configs/best_fit_parameters.toml  final/latest co-fit parameters
#
# The plotting/model implementation still lives in:
#
#   scripts/legacy/plot_yzgo_2d_data_vs_model_legacy.jl
#
# Run from the repo root with:
#
#   julia --project=. scripts/plot_2d_data_vs_model.jl


# ---------------------------------------------------------------------------
# Locate repo root
# ---------------------------------------------------------------------------

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))


# ---------------------------------------------------------------------------
# Load repo utility module
# ---------------------------------------------------------------------------

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit


# ---------------------------------------------------------------------------
# Load configs
# ---------------------------------------------------------------------------

const PLOT_2D_CONTROLS_PATH = joinpath(
    REPO_ROOT,
    "configs",
    "plot_2d_controls.toml",
)

const BEST_FIT_PATH = joinpath(
    REPO_ROOT,
    "configs",
    "best_fit_parameters.toml",
)

const plot_controls = load_toml_config(PLOT_2D_CONTROLS_PATH)
const best_fit_parameters = load_best_fit_parameters(BEST_FIT_PATH)


# ---------------------------------------------------------------------------
# Resolve repo-relative paths
# ---------------------------------------------------------------------------

const NEUTRON_2D_DIR = joinpath(
    REPO_ROOT,
    splitpath(plot_controls["data"]["neutron_2d_subdir"])...,
)

const OUTFIG_DIR = joinpath(
    REPO_ROOT,
    splitpath(plot_controls["output"]["figure_subdir"])...,
)

mkpath(OUTFIG_DIR)


# ---------------------------------------------------------------------------
# Set environment variables expected by the current legacy plotting script
# ---------------------------------------------------------------------------
#
# Important:
#
# The legacy script defines constants using ENV at include time. Therefore
# these ENV values must be set BEFORE including the legacy script.

ENV["YZGO_DATA_DIR"] = NEUTRON_2D_DIR
ENV["YZGO_OUT_DIR"] = OUTFIG_DIR
ENV["MAKIE_BACKEND"] = plot_controls["makie"]["backend"]

# New bridge: make the 2D legacy script read the same best-fit parameter file
# as the co-fit drivers.
ENV["YZGO_BEST_FIT_PARAMETERS_TOML"] = BEST_FIT_PATH


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

println("2D plot controls config:")
println(PLOT_2D_CONTROLS_PATH)
println()

println("Best-fit parameter config:")
println(BEST_FIT_PATH)
println()

println("2D neutron data directory:")
println(NEUTRON_2D_DIR)
println()

println("2D figure output directory:")
println(OUTFIG_DIR)
println()

println("Makie backend:")
println(ENV["MAKIE_BACKEND"])
println()


# ---------------------------------------------------------------------------
# Load legacy 2D plotting implementation
# ---------------------------------------------------------------------------

include(joinpath(
    REPO_ROOT,
    "scripts",
    "legacy",
    "plot_yzgo_2d_data_vs_model_legacy.jl",
))


# ---------------------------------------------------------------------------
# Run 2D comparison
# ---------------------------------------------------------------------------

result_2d = run_yzgo_2d_prelim_model_comparison()


# ---------------------------------------------------------------------------
# Console completion message
# ---------------------------------------------------------------------------

println()
println("2D data-versus-model plot completed.")
println("Controls loaded from:")
println(PLOT_2D_CONTROLS_PATH)
println()
println("Initial/best-fit parameters loaded from:")
println(BEST_FIT_PATH)
println()
println("Wrote figures to:")
println(OUTFIG_DIR)