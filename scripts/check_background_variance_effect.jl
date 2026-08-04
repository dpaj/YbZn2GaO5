#!/usr/bin/env julia
# ACCEPTANCE TEST for the background-variance term: does it stop the fit chasing the 2.08 meV
# magnet artefact?
#
#   julia -t 32 --project=. scripts/check_background_variance_effect.jl
#
# THE FALSIFIABLE CRITERION
#
# The DGX scanned gzz per cut and found the six-cut minimum at 3.60 to be a compromise between one
# artefact and two background-contaminated cuts, outvoting the single physically clean constraint:
#
#   (0,1,0)@9T        prefers 3.40   mode mid-window, exchange cancels at q=0 -- THE CLEAN ONE
#   (0.33,0.33,0)@9T  prefers 3.70   mode at ~1.05 meV, so its 1.8-2.4 meV chi2 is artefact
#   (0.5,0,0)@9T      prefers 3.70   same
#   (0,1,0)@14T       prefers 3.70   ARTEFACT: its mode exits the 3.0 meV window at gzz = 3.702
#
# If the background variance works, the two 9 T DISPERSIVE cuts should move OFF 3.70 toward the
# clean value, because the band driving them there gets large error bars. If they do not move, the
# chain from Stage 0 through Stage 1 to the objective is not doing what it was built to do.
#
# (0,1,0) should be LEAST affected -- three independent lines of evidence say it carries no 2.08 meV
# feature, so its envelope is small and its weights barely change. That is a control, not a null
# result: if (0,1,0) moves a lot, the variance is doing something other than what is intended.
#
# (0,1,0)@14T will NOT be fixed by this. Its problem is the mode leaving the fit window, which is a
# separate defect handled by not using that cut for gzz and by the tightened bounds.
#
# ONE DELIBERATE ASYMMETRY. The intensity scale is fitted with counting errors only in BOTH cases
# (sv_neutron_weighted_scale does not take the background term). That is on purpose: the scale is a
# nuisance parameter, and reweighting it too would confound "the fit stopped chasing the artefact"
# with "the fit renormalised itself".

using Printf, Statistics, LinearAlgebra, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

BLAS.set_num_threads(1)
controls = SV.sv_load_controls(REPO_ROOT)
kc = controls["kpm"]
kc["dims"] = [3, 3, 1]; kc["system_size"] = [36, 36, 1]; kc["repeat_factor"] = [12, 12, 1]
kc["tol"] = 0.05; kc["maxiters"] = 1000; kc["regularization"] = 1e-5
kc["energy_min_meV"] = 0.1; kc["energy_max_meV"] = 4.2; kc["n_energy"] = 241
kc["kernel_fwhm_meV"] = 0.05
kc["thread_max_chunks"] = 32; kc["thread_min_q_per_chunk"] = 1; kc["thread_min_q"] = 1
eh = kc["experimental_histogram"]
eh["enabled"] = true; eh["mode"] = "analytical_cut_volume_grid"
eh["n_measured_h"] = 3; eh["n_measured_k"] = 3; eh["n_measured_l"] = 1
qa = kc["q_averaging"]
qa["enabled"] = true; qa["n_h"] = 3; qa["n_k"] = 3; qa["n_l"] = 1

(; params) = SV.sv_load_params(REPO_ROOT, controls)
params = merge(params, (; J1_meV=0.15, J2_meV=0.01, sigma_J=0.50, sigma_gzz=0.80))
cuts = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
bg = SV.sv_load_background_sigma(REPO_ROOT, cuts)

const NREAL = Int(get(ENV, "YZGO_BGV_NREAL", "") == "" ? 4 : parse(Int, ENV["YZGO_BGV_NREAL"]))
const GZZ = [3.30, 3.40, 3.50, 3.60, 3.70]
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/background_variance_effect")
mkpath(TDIR)

@printf("threads %d, %d realizations, %d cuts, %d q/cut\n", Threads.nthreads(), NREAL,
        length(cuts), length(SV.sv_kpm_1d_q_sampler(cuts[1], controls).qs))
for (k, c) in enumerate(cuts)
    v = filter(x -> x > 0, bg[k])
    @printf("  %-14s %5.1f T  background sigma: median %.3e over %d of %d points\n",
            c.qtag, c.field_T, isempty(v) ? 0.0 : median(v), length(v), length(bg[k]))
end
flush(stdout)

io = open(joinpath(TDIR, "gzz_scan_variance_on_off.csv"), "w")
println(io, "weighting,gzz,chi2_red_total,inflation,field_T,qtag,chi2_red_cut")
flush(io)

results = Dict{Tuple{String,Float64},Any}()
t0 = time()
for g in GZZ, mode in ("off", "on")
    p = merge(params, (; gzz=g))
    o = SV.sv_neutron_objective(p, controls, cuts; realizations=0:(NREAL-1), threaded=true,
            maxiters=1000, relax_attempts=1, on_failure=:record,
            background_sigma = mode == "on" ? bg : nothing)
    results[(mode, g)] = o
    for pc in o.per_cut
        @printf(io, "%s,%.3f,%.6g,%.6g,%.1f,%s,%.6g\n", mode, g, o.chi2_red,
                o.background_variance_inflation, pc.field_T, pc.qtag, pc.chi2_red)
    end
    flush(io)
    @printf("  gzz=%.2f variance %-3s  chi2_red = %9.4f  inflation %.3f  (%.1f h elapsed)\n",
            g, mode, o.chi2_red, o.background_variance_inflation, (time()-t0)/3600)
    flush(stdout)
end
close(io)

println("\n", repeat("=", 78))
println("PER-CUT gzz PREFERENCE, variance OFF vs ON")
println(repeat("=", 78))
tags = unique([(pc.qtag, pc.field_T) for pc in results[("off", GZZ[1])].per_cut])
@printf("%-16s %6s %14s %14s   %s\n", "cut", "field", "prefers OFF", "prefers ON", "verdict")
moved = String[]
for (q, B) in tags
    getc(mode, g) = only(pc.chi2_red for pc in results[(mode, g)].per_cut
                         if pc.qtag == q && pc.field_T ≈ B)
    goff = GZZ[argmin([getc("off", g) for g in GZZ])]
    gon  = GZZ[argmin([getc("on",  g) for g in GZZ])]
    verdict = goff == gon ? "unchanged" : @sprintf("MOVED %+.2f", gon - goff)
    goff != gon && push!(moved, "$q@$(round(Int,B))T")
    @printf("%-16s %5.0fT %14.2f %14.2f   %s\n", q, B, goff, gon, verdict)
end
for mode in ("off", "on")
    g = GZZ[argmin([results[(mode, gg)].chi2_red for gg in GZZ])]
    @printf("\nsix-cut minimum, variance %-3s : gzz = %.2f  (chi2_red %.4f)\n", mode, g,
            results[(mode, g)].chi2_red)
end
println()
println(isempty(moved) ?
  "NO cut changed its preference. The chain is not doing what it was built to do -- check that\n" *
  "the envelope table matches the cuts' energy grids and that bg_sigma is nonzero in 1.8-2.4 meV." :
  "cuts whose preference moved: " * join(moved, ", "))
@printf("\ntotal %.2f h  ->  %s\n", (time()-t0)/3600,
        joinpath(TDIR, "gzz_scan_variance_on_off.csv"))
