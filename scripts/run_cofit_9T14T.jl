# scripts/run_cofit_9T14T.jl
#
# Main YbZn2GaO5 neutron + magnetization co-fit driver.
#
# Run from the repo root with:
#
#     julia --project=. scripts/run_cofit_9T14T.jl

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
    "cofit_9T14T_shared_fraction",
)

result = run_yzgo_neutron_magnetization_cofit_shared_fraction(
    base_dir = NEUTRON_1D_DIR,
    magnetization_csv = MAGNETIZATION_CSV,
    outdir = OUTFIT_DIR,

    run_optimization = false,
    make_plots = true,
    display_figures = true,

    fields_T = [9.0, 14.0],
    neutron_fit_Ei_meV = 4.65,
    neutron_fit_temperature_K = 0.07,
    qtags = ["0_1_0", "0p33_0p33_0", "0p5_0_0"],
    data_mode = :tail_bgsub,

    neutron_weight = 1.0,
    magnetization_weight = 10.0,

    maxiters = 1000,
    n_samples_per_cut = 100_000,
    final_n_samples_per_cut = 500_000,
    n_samples_magnetization = 15_000,
    final_n_samples_magnetization = 150_000,
)

println()
println("Co-fit completed.")
println("Wrote results to:")
println(OUTFIT_DIR)