# scripts/run_cofit_9T14T_smoke.jl
#
# Smoke test for the YbZn2GaO5 neutron + magnetization co-fit.
#
# This intentionally uses:
#   - repo-relative paths,
#   - no optimization,
#   - reduced Monte Carlo sampling,
#   - no interactive figure display.
#
# It reads:
#
#   configs/cofit_controls_smoke.toml   cheap smoke-test controls
#   configs/best_fit_parameters.toml    initial guesses / latest best fit
#
# Run from the repo root with:
#
#   julia --project=. scripts/run_cofit_9T14T_smoke.jl


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
# Load legacy co-fit implementation
# ---------------------------------------------------------------------------
#
# This defines:
#
#   cofit_default_param_specs(...)
#   run_yzgo_neutron_magnetization_cofit_shared_fraction(...)

include(joinpath(
    REPO_ROOT,
    "scripts",
    "legacy",
    "YZGO_cofit_9T14T_shared_fraction_legacy.jl",
))


# ---------------------------------------------------------------------------
# Load configs
# ---------------------------------------------------------------------------

const CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "cofit_controls_smoke.toml")
const BEST_FIT_PATH = joinpath(REPO_ROOT, "configs", "best_fit_parameters.toml")

const controls = load_cofit_controls(CONTROLS_PATH)
const best_fit_parameters = load_best_fit_parameters(BEST_FIT_PATH)
const initial_guess_kwargs = cofit_initial_guess_kwargs(best_fit_parameters)

# Convert TOML initial guesses into legacy co-fit parameter specs.
# This must come AFTER `initial_guess_kwargs` is defined and AFTER the legacy
# script has been included, because `cofit_default_param_specs` lives there.
const specs = cofit_default_param_specs(; initial_guess_kwargs...)


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
# Console summary
# ---------------------------------------------------------------------------

println("Smoke-test controls config:")
println(CONTROLS_PATH)
println()

println("Best-fit parameter config:")
println(BEST_FIT_PATH)
println()

print_initial_guess_kwargs(initial_guess_kwargs)


# ---------------------------------------------------------------------------
# Run smoke test
# ---------------------------------------------------------------------------

result = run_yzgo_neutron_magnetization_cofit_shared_fraction(
    base_dir = NEUTRON_1D_DIR,
    magnetization_csv = MAGNETIZATION_CSV,
    outdir = OUTFIT_DIR,

    # Initial guesses / fit parameter specs
    specs = specs,

    # Smoke-test optimization controls
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
# Console completion message
# ---------------------------------------------------------------------------

println()
println("Smoke test completed.")
println("Controls loaded from:")
println(CONTROLS_PATH)
println()
println("Initial guesses loaded from:")
println(BEST_FIT_PATH)
println()
println("Wrote results to:")
println(OUTFIT_DIR)