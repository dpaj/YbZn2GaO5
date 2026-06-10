# scripts/compare_1d_4p65_3p32_backgrounds.jl
#
# Repo-native driver for the YbZn2GaO5 1D comparison between Ei = 4.65 meV
# and Ei = 3.32 meV background-subtracted scans.
#
# The background/model implementation still lives in:
#
#   scripts/legacy/yzgo_plot_1d_scans_4p65_3p32_compare_legacy.jl
#
# This driver supplies repo-relative paths and plotting controls from:
#
#   configs/background_compare_controls.toml
#
# Run from the repo root with:
#
#   julia --project=. scripts/compare_1d_4p65_3p32_backgrounds.jl


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
# Load controls
# ---------------------------------------------------------------------------

const CONTROLS_PATH = joinpath(
    REPO_ROOT,
    "configs",
    "background_compare_controls.toml",
)

const controls = load_toml_config(CONTROLS_PATH)


# ---------------------------------------------------------------------------
# Resolve repo-relative paths
# ---------------------------------------------------------------------------

const NEUTRON_1D_DIR = joinpath(
    REPO_ROOT,
    splitpath(controls["data"]["neutron_1d_subdir"])...,
)

const OUTFIG_DIR = joinpath(
    REPO_ROOT,
    splitpath(controls["output"]["figure_subdir"])...,
)

mkpath(OUTFIG_DIR)


# ---------------------------------------------------------------------------
# Set environment variables used by the legacy script at include time
# ---------------------------------------------------------------------------

ENV["YZGO_1D_SCAN_DIR"] = NEUTRON_1D_DIR
ENV["MAKIE_BACKEND"] = get(controls["makie"], "backend", "GLMakie")


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

println("1D 4.65/3.32 meV comparison controls config:")
println(CONTROLS_PATH)
println()

println("1D neutron scan directory:")
println(NEUTRON_1D_DIR)
println()

println("Figure output directory:")
println(OUTFIG_DIR)
println()

println("Makie backend:")
println(ENV["MAKIE_BACKEND"])
println()


# ---------------------------------------------------------------------------
# Load legacy implementation
# ---------------------------------------------------------------------------

include(joinpath(
    REPO_ROOT,
    "scripts",
    "legacy",
    "yzgo_plot_1d_scans_4p65_3p32_compare_legacy.jl",
))


# ---------------------------------------------------------------------------
# Run comparison
# ---------------------------------------------------------------------------

bgsub_ylim_values = Float64.(controls["plotting"]["bgsub_ylim"])
length(bgsub_ylim_values) == 2 || error("Expected plotting.bgsub_ylim to have exactly two values")
bgsub_ylim = (bgsub_ylim_values[1], bgsub_ylim_values[2])

if ENV["MAKIE_BACKEND"] != "GLMakie"
    @warn "The legacy 1D comparison script currently uses GLMakie directly; ignoring non-GLMakie backend request" backend=ENV["MAKIE_BACKEND"]
end

GLMakie.activate!()

scans, bgsub_scans, bg_models, fig = main(
    base_dir = NEUTRON_1D_DIR,
    outdir = OUTFIG_DIR,
    save_png = controls["plotting"]["save_png"],
    show_errorbars = controls["plotting"]["show_errorbars"],
    bgsub_ylim = bgsub_ylim,
    make_diagnostic_plot = controls["plotting"]["make_diagnostic_plot"],
    display_figures = controls["plotting"]["display_figures"],
)


# ---------------------------------------------------------------------------
# Console completion message
# ---------------------------------------------------------------------------

println()
println("1D Ei = 4.65 / 3.32 meV background comparison completed.")
println("Controls loaded from:")
println(CONTROLS_PATH)
println()
println("Wrote figures to:")
println(OUTFIG_DIR)
