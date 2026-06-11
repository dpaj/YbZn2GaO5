# scripts/sunny_validate_magnetization_meanfield.jl
#
# Preliminary Sunny.jl validation bridge for the YbZn2GaO5 magnetization model.
#
# This script is intentionally the first Sunny-validation target because it
# checks canonical parameter loading and the two-component magnetization logic
# before the more expensive large-cell/KPM calculations.
#
# Run from the repo root with:
#
#   julia --project=. scripts/sunny_validate_magnetization_meanfield.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

controls = SunnyValidation.sv_load_controls(REPO_ROOT)
result = SunnyValidation.sv_run_meanfield_magnetization(REPO_ROOT; controls)

println()
println("Sunny mean-field magnetization validation completed.")
println("CSV:    ", result.csv_path)
println("Figure: ", result.fig_path)
