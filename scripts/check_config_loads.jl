# Simple smoke test for the repo configuration.
#
# Run from the repo root with:
#
#     julia --project=. scripts/check_config_loads.jl

include(joinpath(@__DIR__, "..", "src", "YZGOCofit.jl"))

using .YZGOCofit

config_path = joinpath(@__DIR__, "..", "configs", "best_fit_parameters.toml")

params = load_best_fit_parameters(config_path)

println("Loaded parameter config:")
println(config_path)
println()

println("Top-level sections:")
for key in keys(params)
    println("  - ", key)
end

println()
println("Full parsed TOML object:")
display(params)