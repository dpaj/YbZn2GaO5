# scripts/check_neutron_cut_loading.jl
#
# Lightweight check that the Sunny validation layer can load the same
# experimental 1D neutron cuts targeted by the analytical co-fit.
#
# Run from the repo root with:
#
#   julia --project=. scripts/check_neutron_cut_loading.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

controls = SunnyValidation.sv_load_controls(REPO_ROOT)
result = SunnyValidation.sv_check_neutron_cut_loading(REPO_ROOT; controls)

println()
println("Neutron 1D cut loading check completed.")
println("Inventory: ", result.inventory_path)
