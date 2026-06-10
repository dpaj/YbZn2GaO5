# scripts/run_cofit_9T14T.jl
#
# Main YbZn2GaO5 neutron + magnetization co-fit driver.
#
# This is the full, non-smoke version of the driver script. It reads
# extrinsic run controls from:
#
#     configs/cofit_controls.toml
#
# The older fitting/model logic is still in:
#
#     scripts/legacy/YZGO_cofit_9T14T_shared_fraction_legacy.jl
#
# Run from the repo root with:
#
#     julia --project=. scripts/run_cofit_9T14T.jl


# ---------------------------------------------------------------------------
# Locate repo root
# ---------------------------------------------------------------------------

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))


# ---------------------------------------------------------------------------
# Load small repo utility module
# ---------------------------------------------------------------------------
#
# This provides:
#   load_cofit_controls(...)
#   toml_symbol(...)
#
# Later, this module will also provide shared model, fitting, plotting,
# and Sunny-validation utilities.

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit


# ---------------------------------------------------------------------------
# Load legacy co-fit implementation
# ---------------------------------------------------------------------------
#
# For now, the real fitting/model code still lives in the legacy script.
# This driver is repo-native: it handles paths and external controls.

include(joinpath(
    REPO_ROOT,
    "scripts",
    "legacy",
    "YZGO_cofit_9T14T_shared_fraction_legacy.jl",
))


# ---------------------------------------------------------------------------
# Load run controls
# ---------------------------------------------------------------------------

const CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "cofit_controls.toml")
const controls = load_cofit_controls(CONTROLS_PATH)


# ---------------------------------------------------------------------------
# Repo-relative input/output paths
# ---------------------------------------------------------------------------

const NEUTRON_1D_DIR = joinpath(
    REPO_ROOT,
    "data",
    "neutron",
    "CNCS_1d_scans",
)

const MAGNETIZATION_CSV = joinpath(
    REPO_ROOT,
    "data",
    "magnetization",
    "YZGO_MvB_black_curve_digitized_visible.csv",
)

const OUTFIT_DIR = joinpath(
    REPO_ROOT,
    "results",
    "fits",
    controls["output"]["fit_name"],
)


# ---------------------------------------------------------------------------
# Run the full co-fit
# ---------------------------------------------------------------------------

result = run_yzgo_neutron_magnetization_cofit_shared_fraction(
    base_dir = NEUTRON_1D_DIR,
    magnetization_csv = MAGNETIZATION_CSV,
    outdir = OUTFIT_DIR,

    # Optimization controls
    run_optimization = controls["optimization"]["run_optimization"],
    maxiters = controls["optimization"]["maxiters"],

    # Plotting controls
    make_plots = controls["plotting"]["make_plots"],
    display_figures = controls["plotting"]["display_figures"],

    # Monte Carlo / disorder sampling controls
    n_samples_per_cut = controls["sampling"]["n_samples_per_cut"],
    final_n_samples_per_cut = controls["sampling"]["final_n_samples_per_cut"],
    n_samples_magnetization = controls["sampling"]["n_samples_magnetization"],
    final_n_samples_magnetization = controls["sampling"]["final_n_samples_magnetization"],

    # Data selection
    fields_T = controls["data"]["fields_T"],
    neutron_fit_Ei_meV = controls["data"]["neutron_fit_Ei_meV"],
    neutron_fit_temperature_K = controls["data"]["neutron_fit_temperature_K"],
    qtags = controls["data"]["qtags"],
    data_mode = toml_symbol(controls["data"]["data_mode"]),

    # Relative weighting of observables
    neutron_weight = controls["weights"]["neutron_weight"],
    magnetization_weight = controls["weights"]["magnetization_weight"],
)


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

println()
println("Co-fit completed.")
println("Controls loaded from:")
println(CONTROLS_PATH)
println()
println("Wrote results to:")
println(OUTFIT_DIR)