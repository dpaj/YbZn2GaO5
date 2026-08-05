#!/usr/bin/env julia
# THE NOISE FLOOR OF THE NEUTRON OBJECTIVE, in the units a fit actually consumes.
#
#   julia -t 32 --project=. scripts/check_realization_scatter_chi2.jl        (~1.6 h)
#
# WHY THIS IS THE LAST MEASUREMENT NEEDED BEFORE A BIG FIT
#
# The realization floor is on record as "12-15%", but that is a SPECTRUM-shape spread. A fit does not
# consume spectra, it consumes `chi2_red`, and nobody has measured the scatter in those units. Without
# it there is no defensible answer to three questions a big fit has to answer:
#
#   * what convergence tolerance is meaningful, rather than chasing noise;
#   * how many realizations to buy, since cost is linear in n and scatter falls as 1/sqrt(n);
#   * whether a claimed improvement between two parameter points is real.
#
# THE SUBTLETY THAT MAKES THIS WORTH DOING PROPERLY. Under common random numbers the realization set
# is FIXED across an optimizer run, so the objective is deterministic in the parameters -- there is no
# per-evaluation noise to average down. What the finite realization set actually does is shift the
# whole `chi2(params)` surface. So the quantity that matters is NOT the scatter in `chi2` at one point.
# It is the scatter in the DIFFERENCE between two points, because CRN correlates the two evaluations
# and the difference is far better determined than either value. That difference scatter is precisely
# the resolution limit on "is point B better than point A", which is the only question an optimizer
# ever asks.
#
# So this measures both, at two nearby parameter points, with the SAME realization sets:
#
#   level scatter       sd of chi2_red(P1) across independent realization sets
#   difference scatter  sd of [chi2_red(P2) - chi2_red(P1)] across the same sets
#   CRN gain            level / difference -- how much the common random numbers buy
#
# Sets are DISJOINT blocks of the realization index, so they are genuinely independent draws:
# n = 4 gives 0:3, 4:7, 8:11, 12:15 and n = 8 gives 0:7, 8:15. Going 4 -> 8 also tests the expected
# 1/sqrt(n) fall, which is what licenses extrapolating to whatever n the big fit ends up using.
#
# The perturbation is deliberately SMALL and along gzz, the direction the optimizer is most sensitive
# to and the one currently pinning against its bound. A perturbation much larger than the noise would
# measure nothing; one much smaller would be all noise. 0.05 in gzz is roughly the step a
# late-stage simplex takes.

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
const BASE = merge(params, (; J1_meV=0.15, J2_meV=0.01, sigma_J=0.50, sigma_gzz=0.80, gzz=3.50))
const PERT = merge(BASE, (; gzz=3.55))
cuts = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)

# Disjoint realization blocks, so the sets are independent draws rather than overlapping samples.
const SETS = [(4, [0:3, 4:7, 8:11, 12:15]), (8, [0:7, 8:15])]
const TDIR = SV.sv_repo_path(REPO_ROOT,
    "results/feature_tables/sunny_validation/realization_scatter_chi2")
mkpath(TDIR)

@printf("threads %d, %d cuts, %d q/cut\n", Threads.nthreads(), length(cuts),
        length(SV.sv_kpm_1d_q_sampler(cuts[1], controls).qs))
@printf("P1: gzz = %.2f    P2: gzz = %.2f (perturbation %.3f)\n",
        BASE.gzz, PERT.gzz, PERT.gzz - BASE.gzz)
println("realization blocks: ", join(string.(vcat([s for (_, s) in SETS]...)), "  "))
flush(stdout)

evaluate(p, rs) = SV.sv_neutron_objective(p, controls, cuts; realizations=rs, threaded=true,
                      maxiters=1000, relax_attempts=1, on_failure=:record)

io = open(joinpath(TDIR, "realization_scatter.csv"), "w")
println(io, "n_real,block,point,gzz,chi2_red_total,field_T,qtag,chi2_red_cut")
flush(io)

results = Dict{Tuple{Int,Int,String},Any}()
t0 = time()
for (n, blocks) in SETS, (bi, rs) in enumerate(blocks), (tag, p) in (("P1", BASE), ("P2", PERT))
    o = evaluate(p, rs)
    results[(n, bi, tag)] = o
    for pc in o.per_cut
        @printf(io, "%d,%d,%s,%.3f,%.6g,%.1f,%s,%.6g\n", n, bi, tag, p.gzz, o.chi2_red,
                pc.field_T, pc.qtag, pc.chi2_red)
    end
    flush(io)
    @printf("  n=%d block %d %s (gzz %.2f): chi2_red = %9.4f   (%.2f h elapsed)\n",
            n, bi, tag, p.gzz, o.chi2_red, (time()-t0)/3600)
    flush(stdout)
end
close(io)

println("\n", repeat("=", 86))
println("THE NUMBERS A FIT CONSUMES")
println(repeat("=", 86))
summary = NamedTuple[]
for (n, blocks) in SETS
    nb = length(blocks)
    c1 = [results[(n, b, "P1")].chi2_red for b in 1:nb]
    c2 = [results[(n, b, "P2")].chi2_red for b in 1:nb]
    d  = c2 .- c1
    lvl = nb > 1 ? std(c1) : NaN
    dif = nb > 1 ? std(d) : NaN
    push!(summary, (; n, nb, mean1 = mean(c1), lvl, meand = mean(d), dif,
                      gain = (isfinite(lvl) && isfinite(dif) && dif > 0) ? lvl / dif : NaN))
    @printf("\nn = %d realizations, %d independent blocks\n", n, nb)
    @printf("  chi2_red(P1) per block   : %s\n", join((@sprintf("%.4f", v) for v in c1), "  "))
    @printf("  chi2_red(P2) per block   : %s\n", join((@sprintf("%.4f", v) for v in c2), "  "))
    @printf("  difference P2 - P1       : %s\n", join((@sprintf("%+.4f", v) for v in d), "  "))
    @printf("  LEVEL scatter   sd(P1)   : %.4f   (%.2f%% of chi2_red = %.3f)\n",
            lvl, 100*lvl/mean(c1), mean(c1))
    @printf("  DIFFERENCE scatter sd(d) : %.4f   <- the resolution limit on 'is P2 better than P1'\n",
            dif)
    @printf("  mean difference          : %+.4f  -> %s\n", mean(d),
            abs(mean(d)) > 2*dif ? "RESOLVED at >2 sd" : "NOT resolved; within 2 sd of zero")
    @printf("  CRN gain (level/diff)    : %.2fx\n", lvl / dif)
end

if length(summary) >= 2
    a, b = summary[1], summary[2]
    println("\n", repeat("-", 86))
    @printf("SCALING with n, %d -> %d (ideal 1/sqrt(n) would give %.3fx):\n",
            a.n, b.n, sqrt(a.n / b.n))
    @printf("  level scatter      %.4f -> %.4f   ratio %.3fx\n", a.lvl, b.lvl, b.lvl / a.lvl)
    @printf("  difference scatter %.4f -> %.4f   ratio %.3fx\n", a.dif, b.dif, b.dif / a.dif)
    # What n would be needed to resolve a target improvement, extrapolating 1/sqrt(n) from n = b.n.
    println("\nRealizations needed to resolve a given chi2_red improvement at 2 sd,")
    println("extrapolating the DIFFERENCE scatter as 1/sqrt(n) from the largest n measured:")
    for target in (0.1, 0.25, 0.5, 1.0, 2.0)
        nreq = b.n * (2 * b.dif / target)^2
        @printf("  resolve %.2f in chi2_red : n >= %6.1f   (~%.1f h per evaluation at %.0f s/realization)\n",
                target, nreq, nreq * 90 / 3600, 90.0)
    end
end

# Per-cut, because the six-cut sum is dominated by whichever cuts fit worst, so a per-cut floor is
# what tells you which cuts can actually discriminate and which are just adding noise.
println("\n", repeat("-", 86))
println("PER-CUT difference scatter at the largest n (which cuts can discriminate?):")
nmax, blocksmax = SETS[end]
nb = length(blocksmax)
tags = [(pc.qtag, pc.field_T) for pc in results[(nmax, 1, "P1")].per_cut]
@printf("  %-14s %6s %12s %12s %12s\n", "cut", "field", "mean chi2", "sd(diff)", "mean diff")
for (q, B) in tags
    get1(b, t) = only(pc.chi2_red for pc in results[(nmax, b, t)].per_cut
                      if pc.qtag == q && pc.field_T ≈ B)
    c1 = [get1(b, "P1") for b in 1:nb]
    d = [get1(b, "P2") - get1(b, "P1") for b in 1:nb]
    @printf("  %-14s %5.0fT %12.3f %12.4f %+12.4f %s\n", q, B, mean(c1),
            nb > 1 ? std(d) : NaN, mean(d),
            nb > 1 && abs(mean(d)) > 2*std(d) ? "  resolved" : "")
end

@printf("\ntotal %.2f h  ->  %s\n", (time()-t0)/3600, joinpath(TDIR, "realization_scatter.csv"))
println("\nCAVEAT: with only 2 blocks at the largest n, sd is a 1-degree-of-freedom estimate and is")
println("itself uncertain by roughly a factor of 2. Treat the n = 4 figures (3 dof) as the reliable")
println("ones and the n = 8 scaling check as indicative.")
