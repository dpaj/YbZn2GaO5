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
# Run from the repo root with:
#
#     julia --project=. scripts/run_cofit_9T14T_smoke.jl


const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(
    REPO_ROOT,
    "scripts",
    "legacy",
    "YZGO_cofit_9T14T_shared_fraction_legacy.jl",
))

const CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "cofit_controls_smoke.toml")
const controls = load_cofit_controls(CONTROLS_PATH)

const NEUTRON_1D_DIR = joinpath(REPO_ROOT, "data", "neutron", "CNCS_1d_scans")

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

result = run_yzgo_neutron_magnetization_cofit_shared_fraction(
    base_dir = NEUTRON_1D_DIR,
    magnetization_csv = MAGNETIZATION_CSV,
    outdir = OUTFIT_DIR,

    run_optimization = controls["optimization"]["run_optimization"],
    make_plots = controls["plotting"]["make_plots"],
    display_figures = controls["plotting"]["display_figures"],

    n_samples_per_cut = controls["sampling"]["n_samples_per_cut"],
    final_n_samples_per_cut = controls["sampling"]["final_n_samples_per_cut"],
    n_samples_magnetization = controls["sampling"]["n_samples_magnetization"],
    final_n_samples_magnetization = controls["sampling"]["final_n_samples_magnetization"],

    fields_T = controls["data"]["fields_T"],
    neutron_fit_Ei_meV = controls["data"]["neutron_fit_Ei_meV"],
    neutron_fit_temperature_K = controls["data"]["neutron_fit_temperature_K"],
    qtags = controls["data"]["qtags"],
    data_mode = toml_symbol(controls["data"]["data_mode"]),

    neutron_weight = controls["weights"]["neutron_weight"],
    magnetization_weight = controls["weights"]["magnetization_weight"],
)

println()
println("Smoke test completed.")
println("Wrote results to:")
println(OUTFIT_DIR)