# scripts/sunny_plot_kpm_2d.jl
#
# Compute a model-only Sunny SpinWaveTheoryKPM intensity map along the configured
# q path and save a 2D PNG + long-form CSV.

repo_root = normpath(joinpath(@__DIR__, ".."))

include(joinpath(repo_root, "src", "YZGOCofit.jl"))
include(joinpath(repo_root, "src", "sunny_validation.jl"))

Main.SunnyValidation.sv_run_kpm_2d_path_map(repo_root)
