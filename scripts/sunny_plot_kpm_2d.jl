# scripts/sunny_plot_kpm_2d.jl
#
# Sunny.jl SpinWaveTheoryKPM 2D comparison against the CNCS 2D experimental cuts.
#
# Style mirrors scripts/plot_2d_data_vs_model.jl / legacy analytical plot:
#   top row    = experimental data
#   bottom row = Sunny KPM model
#   columns    = fields, usually 9 T and 14 T
#
# Run from repo root:
#   julia --project=. scripts/sunny_plot_kpm_2d.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation
using Printf

controls = SunnyValidation.sv_load_controls(REPO_ROOT)
result = SunnyValidation.sv_run_kpm_2d_data_model_comparison(REPO_ROOT; controls)

println()
println("Sunny KPM 2D data/model plot completed.")
println("Figure:")
println(result.fig_path)
