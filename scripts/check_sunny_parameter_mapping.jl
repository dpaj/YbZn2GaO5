# scripts/check_sunny_parameter_mapping.jl
#
# Print and save a diagnostic table showing how the canonical analytical co-fit
# parameters are mapped into the Sunny validation calculators.

repo_root = normpath(joinpath(@__DIR__, ".."))

include(joinpath(repo_root, "src", "YZGOCofit.jl"))
include(joinpath(repo_root, "src", "sunny_validation.jl"))

Main.SunnyValidation.sv_check_sunny_parameter_mapping(repo_root)
