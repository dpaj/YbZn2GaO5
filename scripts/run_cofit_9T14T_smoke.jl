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

include(joinpath(
    REPO_ROOT,
    "scripts",
    "legacy",
    "YZGO_cofit_9T14T_shared_fraction_legacy.jl",
))

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
    "cofit_9T14T_smoke",
)

result = run_yzgo_neutron_magnetization_cofit_shared_fraction(
    base_dir = NEUTRON_1D_DIR,
    magnetization_csv = MAGNETIZATION_CSV,
    outdir = OUTFIT_DIR,

    # Smoke-test mode.
    run_optimization = false,
    make_plots = false,
    display_figures = false,

    # Keep this cheap for the first repo integration test.
    n_samples_per_cut = 5_000,
    final_n_samples_per_cut = 10_000,
    n_samples_magnetization = 2_000,
    final_n_samples_magnetization = 5_000,

    # Same scientific data selection as the co-fit.
    fields_T = [9.0, 14.0],
    neutron_fit_Ei_meV = 4.65,
    neutron_fit_temperature_K = 0.07,
    qtags = ["0_1_0", "0p33_0p33_0", "0p5_0_0"],
    data_mode = :tail_bgsub,
)

println()
println("Smoke test completed.")
println("Wrote results to:")
println(OUTFIT_DIR)