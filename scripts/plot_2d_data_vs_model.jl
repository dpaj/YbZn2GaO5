# scripts/plot_2d_data_vs_model.jl
#
# Repo-native driver for the YbZn2GaO5 2D data-versus-model plot.
#
# The plotting/model implementation still lives in:
#
#   scripts/legacy/plot_yzgo_2d_data_vs_model_legacy.jl
#
# This driver provides clean repo-relative paths and plotting controls.
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
# Load 2D plotting controls
# ---------------------------------------------------------------------------

const PLOT_2D_CONTROLS_PATH = joinpath(
    REPO_ROOT,
    "configs",
    "plot_2d_controls.toml",
)

const plot_controls = load_toml_config(PLOT_2D_CONTROLS_PATH)


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
# The legacy script defines constants such as YZGO_2D_DATA_DIR using ENV at
# include time. Therefore these ENV values must be set BEFORE including the
# legacy script.

ENV["YZGO_DATA_DIR"] = NEUTRON_2D_DIR
ENV["YZGO_OUT_DIR"] = OUTFIG_DIR
ENV["MAKIE_BACKEND"] = plot_controls["makie"]["backend"]


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

println("2D plot controls config:")
println(PLOT_2D_CONTROLS_PATH)
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
#
# This include must happen after setting ENV["YZGO_DATA_DIR"],
# ENV["YZGO_OUT_DIR"], and ENV["MAKIE_BACKEND"].

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
println("Wrote figures to:")
println(OUTFIG_DIR)