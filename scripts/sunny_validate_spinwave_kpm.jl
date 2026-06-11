# scripts/sunny_validate_spinwave_kpm.jl
#
# Preliminary Sunny.jl SpinWaveTheoryKPM validation for the YbZn2GaO5 1D cuts.
#
# This is a first real scaffold for the eventual neutron validation:
#
#   total = Sunny dispersive component + r2 * Sunny flat component
#
# It uses the same qtags/fields as the analytical co-fit target.  The intensity
# extraction from Sunny's result object may need minor API tuning depending on
# the installed Sunny version; if it errors, inspect the printed propertynames.
#
# Run from the repo root with:
#
#   julia --project=. scripts/sunny_validate_spinwave_kpm.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

controls = SunnyValidation.sv_load_controls(REPO_ROOT)
result = SunnyValidation.sv_run_kpm_1d(REPO_ROOT; controls)

println()
println("Sunny KPM 1D validation completed.")
println("Figure: ", result.fig_path)
