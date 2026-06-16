# scripts/extract_neutron_magnetization_features.jl
#
# Empirical feature extraction for YbZn2GaO5 neutron 1D scans and magnetization.
#
# Run from the repo root with:
#
#     julia --project=. scripts/extract_neutron_magnetization_features.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit
using CairoMakie

const CONTROLS_PATH = joinpath(REPO_ROOT, "configs", "feature_extraction_controls.toml")
const controls = load_toml_config(CONTROLS_PATH)

println("Feature extraction controls config:")
println(CONTROLS_PATH)
println()

# This first modern feature-extraction workflow uses CairoMakie for static PNG
# output. Keep this non-interactive for reproducibility and easy CI-style checks.
CairoMakie.activate!()

result = run_neutron_magnetization_feature_extraction(
    repo_root = REPO_ROOT,
    controls = controls,
)
