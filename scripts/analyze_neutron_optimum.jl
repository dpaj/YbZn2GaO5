#!/usr/bin/env julia
# Follow-up analysis once the multi-start Nelder-Mead runs have finished.
#
#   julia -t auto --project=. scripts/analyze_neutron_optimum.jl
#
# Five stages, ordered so that the robust ones come last and therefore still run if an
# earlier one overruns its deadline:
#
#   1  post-cache intra-process thread curve -- the existing numbers predate the KPM
#      operator cache and are stale; this also sets the launch recipe from here on
#   2  collect every start, rank them, and report whether they agree (multiple minima?)
#   3  profile likelihood in each free parameter -> ACTUAL ERROR BARS
#   4  seed validation -- is the optimum real, or an artefact of one disorder draw?
#   5  M(H) at the neutron optimum -- cross-observable check with NO choice of weights
#
# Stage 5 is the scientifically interesting one: M(H) constrains B_sat ~ J1/gzz, and the
# fitted parameters move that ratio down by roughly a third from the by-eye pair, i.e.
# toward what M(H) independently prefers. If M(H) improves at parameters fitted only to
# neutrons, that is cross-observable evidence obtained without weighting anything.

using Printf, Statistics, LinearAlgebra, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const LOADED = SV.sv_load_diagnostic_controls(REPO_ROOT,
    "configs/neutron_optimization_controls.toml"; env_var="YZGO_NEUTRON_OPT_CONTROLS")
const CFG = LOADED.diag
const RUN = get(CFG, "run", Dict{String,Any}())
const FU = get(CFG, "followup", Dict{String,Any}())
const controls = LOADED.controls

params0, _ = SV.sv_apply_param_overrides(
    SV.sv_load_params(REPO_ROOT, controls).params, RUN)
const CUTS = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
const NREAL = Int(get(RUN, "n_realizations", 4))
const MAXITERS = Int(get(RUN, "minimize_maxiters", 1000))
const RELAX = Int(get(RUN, "relax_attempts", 1))
const FREE = Symbol.(String.(get(RUN, "free_parameters",
    ["J1_meV", "sigma_J", "gzz", "sigma_gzz"])))
const LO = Float64.(get(RUN, "lower_bounds", [0.05, 0.0, 2.8, 0.0]))
const HI = Float64.(get(RUN, "upper_bounds", [0.40, 0.60, 4.2, 1.4]))
const TDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "table_subdir",
    "results/feature_tables/sunny_validation/neutron_optimization"))
mkpath(TDIR)

const T0 = time()
wall_h() = (time() - T0) / 3600
setp(p, xs) = merge(p, NamedTuple{Tuple(FREE)}(Tuple(xs)))

banner(s) = (@printf("\n%s\n%s\n%s\n", repeat("=", 78), s, repeat("=", 78)); flush(stdout))

function neutron(p; nreal=NREAL, cuts=CUTS)
    try
        return SV.sv_neutron_objective(p, controls, cuts; realizations=0:(nreal - 1),
            threaded=true, maxiters=MAXITERS, relax_attempts=RELAX, on_failure=:record)
    catch err
        @printf("    !! failed: %s\n", first(split(sprint(showerror, err), '\n'))[1:min(60, end)])
        flush(stdout)
        return nothing
    end
end

row!(io, xs...) = (println(io, join(xs, ",")); flush(io))

# ---------------------------------------------------------------- 1  thread curve
banner("Stage 1 -- post-cache intra-process thread curve")
println("""
The chunk-scaling numbers in CLAUDE.md predate the KPM operator cache (7f8016b). The DGX
saw their peak move 5.87x -> 13.51x with the knee shifting 16 -> 64 threads once
construction was hoisted out of the threaded loop, so this box's curve is stale and the
launch recipe may have changed. Measured on ONE cut so it is cheap.
""")
open(joinpath(TDIR, "thread_curve_postcache.csv"), "w") do io
    row!(io, "chunks", "kpm_seconds", "speedup", "efficiency_pct")
    gcut = filter(c -> c.qtag == "0_1_0" && c.field_T ≈ 9.0, CUTS)
    neutron(params0; nreal=1, cuts=gcut)          # warm up JIT
    base = 0.0
    for nch in Int.(get(FU, "thread_curve_chunks", [1, 2, 4, 8, 16, 32]))
        controls["kpm"]["thread_max_chunks"] = nch
        o = neutron(params0; nreal=1, cuts=gcut)
        o === nothing && continue
        s = o.kpm_seconds
        nch == 1 && (base = s)
        sp = base / s
        @printf("  %3d chunks  %7.1f s  %6.2fx  %4.0f%%\n", nch, s, sp, 100 * sp / nch)
        row!(io, nch, round(s; digits=2), round(sp; digits=3), round(100 * sp / nch; digits=1))
    end
end
delete!(controls["kpm"], "thread_max_chunks")     # restore default chunking

# ---------------------------------------------------------------- 2  collect starts
banner("Stage 2 -- collect the multi-start results")
bests = NamedTuple[]
for f in sort(readdir(TDIR))
    (startswith(f, "best_") && endswith(f, ".csv")) || continue
    lines = filter(!isempty, strip.(readlines(joinpath(TDIR, f))))
    length(lines) >= 2 || continue
    hdr = String.(split(lines[1], ',')); val = String.(split(lines[2], ','))
    d = Dict(zip(hdr, val))
    g(k) = something(tryparse(Float64, get(d, k, "")), NaN)
    push!(bests, (; tag = get(d, "tag", f), chi2 = g("chi2_red"),
                    x = [g(String(s)) for s in FREE], nev = g("n_evals"),
                    hours = g("hours"), stopped = get(d, "stopped_by", "?")))
end
if isempty(bests)
    error("No best_*.csv in $TDIR -- did the Nelder-Mead starts run?")
end
sort!(bests; by = b -> b.chi2)
@printf("%-14s %12s %8s %7s  %s\n", "start", "chi2_red", "evals", "hours", "parameters")
for b in bests
    @printf("%-14s %12.5g %8.0f %7.2f  %s\n", b.tag, b.chi2, b.nev, b.hours,
            join([@sprintf("%s=%.4f", s, v) for (s, v) in zip(FREE, b.x)], " "))
end
open(joinpath(TDIR, "multistart_summary.csv"), "w") do io
    row!(io, "rank", "tag", join(String.(FREE), ","), "chi2_red", "n_evals", "hours", "stopped_by")
    for (i, b) in enumerate(bests)
        row!(io, i, b.tag, join(round.(b.x; digits=5), ","), round(b.chi2; sigdigits=8),
             Int(b.nev), round(b.hours; digits=3), b.stopped)
    end
end
xbest = bests[1].x
pbest = setp(params0, xbest)
@printf("\nbest: %s  chi2_red = %.5g\n",
        join([@sprintf("%s=%.4f", s, v) for (s, v) in zip(FREE, xbest)], " "), bests[1].chi2)
# Do the starts agree? Spread across the top few is the real test for multiple minima.
top = [b for b in bests if isfinite(b.chi2) && b.chi2 < 1.5 * bests[1].chi2]
if length(top) > 1
    println("\nspread across starts within 1.5x of the best (multiple-minimum check):")
    for (j, s) in enumerate(FREE)
        v = [b.x[j] for b in top]
        @printf("  %-11s %.4f - %.4f  (spread %.1f%% of the best value)\n", String(s),
                minimum(v), maximum(v), 100 * (maximum(v) - minimum(v)) / max(1e-9, abs(xbest[j])))
    end
    println(maximum(100 * (maximum(b.x[j] for b in top) - minimum(b.x[j] for b in top)) /
                    max(1e-9, abs(xbest[j])) for j in eachindex(FREE)) < 10 ?
        "  => starts AGREE; the optimum looks unique in this box." :
        "  => starts DISAGREE; there are multiple minima and the reported optimum is one of several.")
end

# ---------------------------------------------------------------- 3  confirm at 225 q
banner("Stage 3 -- confirm the optimum at the properly converged q setting")
ms = Int(get(FU, "confirm_measured_side", 5)); rs = Int(get(FU, "confirm_resolution_side", 3))
o81 = neutron(pbest)
eh = controls["kpm"]["experimental_histogram"]; qa = controls["kpm"]["q_averaging"]
eh["n_measured_h"] = ms; eh["n_measured_k"] = ms
qa["n_h"] = rs; qa["n_k"] = rs
o225 = neutron(pbest)
if o81 !== nothing && o225 !== nothing
    @printf("  81 q  (3x3 x 3x3): chi2_red = %.6g\n", o81.chi2_red)
    @printf("  %d q (%dx%d x %dx%d): chi2_red = %.6g   (%.2f%% change)\n",
            ms^2 * rs^2, ms, ms, rs, rs, o225.chi2_red,
            100 * (o225.chi2_red - o81.chi2_red) / o81.chi2_red)
    println("  The DGX measured 81 q as 0.64% from converged and 5x5 x 3x3 as 0.058%.")
    println("  Both sit far below the 12-15% realization floor and the background systematic.")
end
eh["n_measured_h"] = 3; eh["n_measured_k"] = 3; qa["n_h"] = 3; qa["n_k"] = 3

# ---------------------------------------------------------------- 4  profile likelihood
banner("Stage 4 -- profile likelihood: error bars for each free parameter")
println("""
Each parameter is scanned with the others HELD at the optimum. That is a conditional
profile, not a full re-optimized one, so the widths it gives are LOWER BOUNDS on the true
uncertainty -- correlations between parameters can only widen them. Judge significance
against the realization floor, not against the curvature alone.
""")
npts = Int(get(FU, "profile_points", 7))
frac = Float64(get(FU, "profile_span_fraction", 0.45))
open(joinpath(TDIR, "profile_likelihood.csv"), "w") do io
    row!(io, "parameter", "value", "chi2_red", "delta_chi2", "rms", "seconds")
    for (j, s) in enumerate(FREE)
        c = xbest[j]
        span = max(frac * abs(c), 0.05 * (HI[j] - LO[j]))
        vals = unique(clamp.(collect(range(c - span, c + span; length = npts)), LO[j], HI[j]))
        @printf("  %s (optimum %.4f):\n", String(s), c)
        for v in vals
            wall_h() > Float64(get(FU, "wall_clock_budget_hours", 20.0)) &&
                (println("    !! budget reached; truncating profiles."); break)
            x = copy(xbest); x[j] = v
            t = time(); o = neutron(setp(params0, x))
            c2 = o === nothing ? NaN : o.chi2_red
            @printf("    %-11s = %8.4f  chi2 = %10.5g  (dchi2 = %+8.4g)\n",
                    String(s), v, c2, c2 - bests[1].chi2)
            row!(io, String(s), round(v; digits=5), round(c2; sigdigits=8),
                 round(c2 - bests[1].chi2; sigdigits=6),
                 o === nothing ? NaN : round(o.rms; sigdigits=6), round(time() - t; digits=1))
        end
    end
end

# ---------------------------------------------------------------- 5  seed validation
banner("Stage 5 -- seed validation: is the optimum real or one disorder draw?")
println("""
Every neutron result so far uses seed 20260611 with realizations 0-3. The M(H) protocol
requires validating an optimum against different seeds; the neutron side never has. If
chi2 moves by more than the realization floor when only the disorder DRAW changes, the
optimum is a common-random-numbers artefact rather than physics.
""")
seeds = Int.(get(FU, "validation_seeds", [20260611]))
vnr = Int(get(FU, "validation_realizations", 8))
open(joinpath(TDIR, "seed_validation.csv"), "w") do io
    row!(io, "seed", "n_realizations", "chi2_red_optimum", "chi2_red_byeye", "seconds")
    pbye = merge(params0, (; J1_meV=0.25, sigma_J=0.50, gzz=3.80, sigma_gzz=0.80))
    for sd in seeds
        wall_h() > Float64(get(FU, "wall_clock_budget_hours", 20.0)) && break
        controls["common"]["seed"] = sd
        t = time()
        oa = neutron(pbest; nreal=vnr)
        ob = neutron(pbye; nreal=vnr)
        ca = oa === nothing ? NaN : oa.chi2_red
        cb = ob === nothing ? NaN : ob.chi2_red
        @printf("  seed %-10d optimum chi2 = %10.5g   by-eye chi2 = %10.5g\n", sd, ca, cb)
        row!(io, sd, vnr, round(ca; sigdigits=8), round(cb; sigdigits=8),
             round(time() - t; digits=1))
    end
end
controls["common"]["seed"] = 20260611

# ---------------------------------------------------------------- 6  M(H) cross-check
banner("Stage 6 -- M(H) at the neutron optimum, with no weights chosen")
println("""
M(H) constrains B_sat ~ J1/gzz. The by-eye pair gives ~5.1 T where M(H) independently wants
~4.0 T, and the fitted pair lowers that ratio by roughly a third. So M(H) should IMPROVE at
parameters fitted only to neutrons. That is a genuine cross-observable test, and it needs no
relative weighting -- which is the part of a co-fit that cannot be justified by assertion.
""")
Bs = collect(range(Float64(get(FU, "mvh_field_min_T", 0.2)),
                   Float64(get(FU, "mvh_field_max_T", 6.8));
                   length = Int(get(FU, "mvh_n_fields", 24))))
M_exp = SV.sv_mvh_target(REPO_ROOT, controls, Bs)
mcell = Tuple(Int.(get(FU, "mvh_cell_size", [12, 12, 1])))
mnr = Int(get(FU, "mvh_realizations", 8))
open(joinpath(TDIR, "mvh_cross_check.csv"), "w") do io
    row!(io, "set,J1_meV,sigma_J,gzz,sigma_gzz,J1_over_gzz,rms_uB,max_abs_uB,A_M,chi_vv,nan_fraction")
    for (nm, p) in (("by-eye", merge(params0, (; J1_meV=0.25, sigma_J=0.50, gzz=3.80,
                                                 sigma_gzz=0.80))),
                    ("neutron-optimum", pbest))
        try
            o = SV.sv_mvh_objective(p, controls, Bs, M_exp; cell_size=mcell,
                    seed_dims=(3, 3, 1), realizations=0:(mnr - 1), maxiters=2000,
                    threaded=true)
            @printf("  %-16s rms = %.5f uB   max|res| = %.5f   A_M = %.4g  chi_vv = %.4g\n",
                    nm, o.rms, o.max_abs, o.A_M, o.chi_vv)
            row!(io, nm, round(p.J1_meV; digits=4), round(p.sigma_J; digits=4),
                 round(p.gzz; digits=4), round(p.sigma_gzz; digits=4),
                 round(p.J1_meV / p.gzz; sigdigits=5), round(o.rms; sigdigits=6),
                 round(o.max_abs; sigdigits=6), round(o.A_M; sigdigits=6),
                 round(o.chi_vv; sigdigits=6), round(o.nan_fraction; digits=4))
        catch err
            @printf("  %-16s FAILED: %s\n", nm,
                    first(split(sprint(showerror, err), '\n'))[1:min(70, end)])
        end
        flush(stdout)
    end
end

@printf("\ntotal wall clock %.2f h\ntables -> %s\n", wall_h(), TDIR)
println("done")
