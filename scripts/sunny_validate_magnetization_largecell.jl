# scripts/sunny_validate_magnetization_largecell.jl
#
# Preliminary Sunny.jl large-cell minimize_energy! magnetization validation.
#
# This uses the same canonical best-fit parameters as the analytical co-fit and
# constructs two Sunny components:
#
#   total = magnetization_global_scale * [(Sunny dispersive + r2*Sunny flat)/(1+r2) + chi_vv B]
#
# The first version is deliberately small/cheap.  Increase dims/repeat_factor in
# configs/sunny_validation_controls.toml only after the basic calculation works.
#
# Run from the repo root with:
#
#   julia --project=. scripts/sunny_validate_magnetization_largecell.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

controls = SunnyValidation.sv_load_controls(REPO_ROOT)
result = SunnyValidation.sv_run_largecell_magnetization(REPO_ROOT; controls)

println()
println("Sunny large-cell magnetization validation completed.")
println("CSV:    ", result.csv_path)
println("Figure: ", result.fig_path)
