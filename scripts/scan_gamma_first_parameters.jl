#!/usr/bin/env julia
# Factorized "Gamma-first" parameter scan for the 1D neutron cuts.
#
# THE PHYSICS THIS EXPLOITS
# -------------------------
# At the zone center the exchange cancels identically:
#
#     omega(q) = gzz*mu_B*B + S*[J(0) - J(q)]     ->    omega(0) = gzz*mu_B*B
#
# so the (0,1,0) cut position depends on gzz ALONE and its width on sigma_gzz alone.
# g-factor disorder broadens everywhere in the zone; exchange disorder only "turns on"
# away from q = 0. That gives a factorization:
#
#   Stage 1   (0,1,0) at 9 and 14 T          ->  gzz (peak position), sigma_gzz (width)
#   Stage 2   K and M at 9 and 14 T          ->  J1 (bandwidth), sigma_J (excess width)
#   Stage 3   (0,1,0) again at stage-2 J     ->  did stage 1 move? (consistency)
#
# WHY THIS IS NOT AN IDENTITY, AND WHAT STAGE 1B MEASURES
# ------------------------------------------------------
# The measured "Gamma" cut is a VOLUME, not a point: n_measured x n_resolution q
# samples spread around (0,1,0). J therefore leaks in through the finite q extent.
# Stage 1b evaluates the stage-1 optimum across a range of J1 to quantify that leak,
# so the factorization is validated rather than assumed. If stage 1b shows the Gamma
# chi2 moving appreciably with J1, the factorization is only a starting point and the
# stages must be iterated.
#
# THE HYPOTHESIS UNDER TEST IN STAGE 2
# ------------------------------------
# An earlier sigma_J scan at FIXED J1 showed chi2 falling monotonically with sigma_J.
# The suspicion is that this is an artifact: the model bandwidth at K and M is 10-20%
# too small, and broadening is the only handle available to push spectral weight up to
# where the data want it. If so, the stage-2 chi2 valley will run diagonally from
# (low J1, high sigma_J) to (higher J1, low sigma_J), and the true optimum will sit at
# larger J1 and modest sigma_J. A diagonal valley is the signature of sigma_J acting as
# a bandwidth surrogate; a valley aligned with the sigma_J axis would mean the
# broadening is real and independent.
#
# UNATTENDED-RUN DESIGN
# ---------------------
# This is built to run overnight with no supervision:
#   * every row is appended to CSV the moment it is computed, so a kill at any point
#     leaves usable partial results;
#   * a calibration point is timed FIRST and the projected total printed, so the log
#     says up front whether the job fits in the budget;
#   * a wall-clock budget skips later stages rather than overrunning;
#   * a failed parameter point is recorded and skipped, never fatal.
#
# Run with threads (KPM threads over q):
#   julia -t auto --project=. scripts/scan_gamma_first_parameters.jl

using Printf, Statistics, Dates, LinearAlgebra, Sunny

const REPO_ROOT = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

# ---------------------------------------------------------------- controls

const LOADED = SV.sv_load_diagnostic_controls(REPO_ROOT,
    "configs/gamma_first_scan_controls.toml"; env_var="YZGO_GAMMA_FIRST_CONTROLS")
const CFG = LOADED.diag
const RUN = get(CFG, "run", Dict{String,Any}())
const controls = LOADED.controls   # control_overrides are already deep-merged in

# sv_apply_param_overrides returns a TUPLE; destructuring it is not optional.
params, applied = SV.sv_apply_param_overrides(
    SV.sv_load_params(REPO_ROOT, controls).params, RUN)

const CUTS_ALL = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
const GAMMA_TAGS = String.(get(RUN, "gamma_qtags", ["0_1_0"]))
const DISP_TAGS = String.(get(RUN, "dispersive_qtags", ["0p33_0p33_0", "0p5_0_0"]))
cuts_gamma = filter(c -> c.qtag in GAMMA_TAGS, CUTS_ALL)
cuts_disp = filter(c -> c.qtag in DISP_TAGS, CUTS_ALL)

const NREAL = Int(get(RUN, "n_realizations", 8))
const REALS = 0:(NREAL - 1)
const MAXITERS = Int(get(RUN, "minimize_maxiters", 1000))
const BUDGET_H = Float64(get(RUN, "wall_clock_budget_hours", 11.0))
const OUTDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "table_subdir",
    "results/feature_tables/sunny_validation/gamma_first_scan"))
mkpath(OUTDIR)

grid(key, default) = Float64.(get(RUN, key, default))
const GZZ_GRID = grid("gzz_grid", [3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
const SGZZ_GRID = grid("sigma_gzz_grid", [0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2])
const J1_GRID = grid("J1_grid", [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45])
const SJ_GRID = grid("sigma_J_grid", [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
const J1_PROBE = grid("J1_probe_grid", [0.15, 0.20, 0.25, 0.30, 0.35, 0.40])

const T0 = time()
elapsed_h() = (time() - T0) / 3600
remaining_h() = BUDGET_H - elapsed_h()

# Per-stage deadlines as a fraction of the budget. Without these, a stage-1 grid that
# runs slower than projected silently consumes the whole night and stage 4 -- the stage
# that actually protects against the sequential-factorization trap -- never runs at all.
# Truncating stage 1 costs a few grid points, all of whose rows are already on disk;
# losing stage 4 costs the only guard against reporting a globally worse optimum.
const DEADLINES = Dict(
    "stage1"  => 0.30, "stage1b" => 0.35, "stage2" => 0.65,
    "stage3"  => 0.72, "stage4"  => 0.95,
)
past_deadline(stage) = elapsed_h() > DEADLINES[stage] * BUDGET_H
deadline_h(stage) = DEADLINES[stage] * BUDGET_H

function banner(s)
    @printf("\n%s\n%s\n%s\n", repeat("=", 78), s, repeat("=", 78))
    flush(stdout)
end

# ---------------------------------------------------------------- evaluation

# One parameter point against one cut subset. Errors are recorded, never fatal: an
# unattended 12 h job must not die on a single KPM abort.
function evaluate(p, cuts; label="")
    t = time()
    try
        o = SV.sv_neutron_objective(p, controls, cuts;
            realizations=REALS, threaded=true, maxiters=MAXITERS, on_failure=:record)
        return (; ok=o.ok, chi2=o.chi2_red, rms=o.rms, scale=o.scale,
                  n_failed=o.n_failed, seconds=time() - t,
                  reg=isempty(o.regularization_values) ? NaN : maximum(o.regularization_values),
                  per_cut=o.per_cut)
    catch err
        msg = first(split(sprint(showerror, err), '\n'))
        @printf("    !! %s failed: %s\n", label, msg[1:min(70, end)])
        flush(stdout)
        return (; ok=false, chi2=Inf, rms=Inf, scale=NaN, n_failed=-1,
                  seconds=time() - t, reg=NaN, per_cut=NamedTuple[])
    end
end

# Append-as-you-go CSV. Written per row so a kill leaves usable partial output.
function open_csv(name, header)
    path = joinpath(OUTDIR, name)
    io = open(path, "w")
    println(io, header)
    flush(io)
    return (io, path)
end

function row!(io, fields...)
    println(io, join(fields, ","))
    flush(io)
end

# ---------------------------------------------------------------- calibration

banner("Gamma-first factorized parameter scan")
@printf("threads              %d\n", Threads.nthreads())
@printf("realizations         %d  (common random numbers, realization = 0:%d)\n",
        NREAL, NREAL - 1)
@printf("cell                 %s\n", string(get(controls["kpm"], "system_size", "?")))
@printf("kpm tol              %s\n", string(get(controls["kpm"], "tol", "?")))
@printf("regularization       %s (pinned)\n", string(get(controls["kpm"], "regularization", "?")))
@printf("minimize maxiters    %d\n", MAXITERS)
@printf("Gamma cuts           %d  %s\n", length(cuts_gamma),
        join(["$(c.qtag)@$(c.field_T)T" for c in cuts_gamma], " "))
@printf("dispersive cuts      %d  %s\n", length(cuts_disp),
        join(["$(c.qtag)@$(c.field_T)T" for c in cuts_disp], " "))
@printf("q per cut            %d\n", length(SV.sv_kpm_1d_q_sampler(CUTS_ALL[1], controls).qs))
@printf("wall-clock budget    %.1f h\n", BUDGET_H)
isempty(applied) || @printf("param overrides      %s\n", join(applied, ", "))

banner("Calibration -- timing one parameter point before committing to the grids")
cal_g = evaluate(params, cuts_gamma; label="calib-gamma")
cal_d = evaluate(params, cuts_disp; label="calib-disp")
@printf("  Gamma subset      %6.1f s/point   chi2_red = %-10.4g scale = %.4g\n",
        cal_g.seconds, cal_g.chi2, cal_g.scale)
@printf("  dispersive subset %6.1f s/point   chi2_red = %-10.4g scale = %.4g\n",
        cal_d.seconds, cal_d.chi2, cal_d.scale)

n1 = length(GZZ_GRID) * length(SGZZ_GRID)
n1b = length(J1_PROBE)
n2 = length(J1_GRID) * length(SJ_GRID)
n3 = 25
n4 = 1 + 8 * Int(get(RUN, "refine_sweeps", 2))       # start + 4 params x 2 signs x sweeps
cal_all = cal_g.seconds + cal_d.seconds              # all six cuts, less ground-state reuse
proj = (n1 * cal_g.seconds + n1b * cal_g.seconds + n2 * cal_d.seconds +
        n3 * cal_g.seconds + (n4 + 3) * cal_all) / 3600
@printf("\n  projected total: stage1 %d + stage1b %d + stage2 %d + stage3 %d + stage4 %d = %.1f h\n",
        n1, n1b, n2, n3, n4, proj)
@printf("  budget %.1f h -- %s\n", BUDGET_H,
        proj <= BUDGET_H ? "FITS" : "OVER BUDGET, later stages will be skipped")
flush(stdout)

# ---------------------------------------------------------------- stage 1

banner("Stage 1 -- (0,1,0) only: gzz from peak POSITION, sigma_gzz from WIDTH")
println("At q = 0 the exchange cancels, so this subset should be insensitive to J1.")
println("The fitted intensity scale is free here, which is what removes amplitude from")
println("the position/width determination -- position and width are scale-invariant.\n")

io1, p1 = open_csv("stage1_gamma_gzz_sigma_gzz.csv",
    "gzz,sigma_gzz,chi2_red,rms,scale,ok,n_failed,seconds,regularization")
best1 = (; chi2=Inf, gzz=params.gzz, sigma_gzz=params.sigma_gzz)
done = 0
for gz in GZZ_GRID, sg in SGZZ_GRID
    global best1, done
    if past_deadline("stage1")
        @printf("  !! stage-1 deadline (%.1f h) reached after %d/%d points; truncating so
",
                deadline_h("stage1"), done, n1)
        println("     the later stages still run. Completed rows are already on disk.")
        break
    end
    p = merge(params, (; gzz=gz, sigma_gzz=sg))
    r = evaluate(p, cuts_gamma; label=@sprintf("gzz=%.2f sgzz=%.2f", gz, sg))
    row!(io1, gz, sg, r.chi2, r.rms, r.scale, r.ok, r.n_failed,
         round(r.seconds; digits=2), r.reg)
    if r.ok && r.chi2 < best1.chi2
        best1 = (; chi2=r.chi2, gzz=gz, sigma_gzz=sg)
    end
    done += 1
    @printf("  [%3d/%3d] gzz=%.2f sigma_gzz=%.2f  chi2=%-11.4g %s (%.0f s, %.1f h left)\n",
            done, n1, gz, sg, r.chi2, r.ok ? " " : "FAIL", r.seconds, remaining_h())
    flush(stdout)
end
close(io1)
@printf("\n  stage 1 best: gzz = %.2f, sigma_gzz = %.2f, chi2_red = %.4g\n",
        best1.gzz, best1.sigma_gzz, best1.chi2)
@printf("  Zeeman check: gzz*mu_B*B = %.3f meV at 9 T, %.3f meV at 14 T\n",
        best1.gzz * 0.05788 * 9, best1.gzz * 0.05788 * 14)
println("  -> $p1")

# ---------------------------------------------------------------- stage 1b

banner("Stage 1b -- is the factorization actually valid? Gamma chi2 vs J1")
println("The Gamma CUT is a q volume, not a point, so J leaks in through its extent.")
println("A flat row here means the factorization holds and stage 2 can be trusted;")
println("appreciable curvature means the stages must be iterated.\n")

io1b, p1b = open_csv("stage1b_gamma_vs_J1.csv",
    "J1_meV,gzz,sigma_gzz,chi2_red,rms,scale,ok,seconds")
probe = NamedTuple[]
for j1 in J1_PROBE
    if past_deadline("stage1b")
        println("  !! stage-1b deadline reached; truncating the J1 probe.")
        break
    end
    p = merge(params, (; gzz=best1.gzz, sigma_gzz=best1.sigma_gzz, J1_meV=j1))
    r = evaluate(p, cuts_gamma; label=@sprintf("J1=%.3f", j1))
    row!(io1b, j1, best1.gzz, best1.sigma_gzz, r.chi2, r.rms, r.scale, r.ok,
         round(r.seconds; digits=2))
    push!(probe, (; j1, chi2=r.chi2))
    @printf("  J1 = %.3f meV   chi2_red = %-11.4g (%.0f s)\n", j1, r.chi2, r.seconds)
    flush(stdout)
end
close(io1b)
if length(probe) > 1
    fin = [x for x in probe if isfinite(x.chi2)]
    if length(fin) > 1
        lo, hi = minimum(x.chi2 for x in fin), maximum(x.chi2 for x in fin)
        @printf("\n  Gamma chi2 varies by %.1f%% across J1 = %.2f-%.2f meV\n",
                100 * (hi - lo) / lo, minimum(x.j1 for x in fin), maximum(x.j1 for x in fin))
        println((hi - lo) / lo < 0.05 ?
            "  => FACTORIZATION HOLDS. Gamma is a clean gzz/sigma_gzz probe; stage 2 is safe." :
            "  => FACTORIZATION LEAKS. The Gamma cut volume carries real J sensitivity,\n" *
            "     so treat stage 1 as a starting point and iterate stages 1 and 2.")
    end
end
println("  -> $p1b")

# ---------------------------------------------------------------- stage 2

banner("Stage 2 -- K and M only: J1 from BANDWIDTH, sigma_J from EXCESS width")
println("gzz and sigma_gzz are now FIXED at the stage-1 values, so the broadening that")
println("g disorder already accounts for is not re-fit here. A diagonal chi2 valley")
println("means sigma_J was standing in for a too-small bandwidth.\n")

best2 = (; chi2=Inf, J1=params.J1_meV, sigma_J=params.sigma_J)
if past_deadline("stage2")
    println("  !! past the stage-2 deadline already; SKIPPING stage 2.")
else
    io2, p2 = open_csv("stage2_disp_J1_sigma_J.csv",
        "J1_meV,sigma_J,gzz,sigma_gzz,chi2_red,rms,scale,ok,n_failed,seconds,regularization")
    done2 = 0
    for j1 in J1_GRID, sj in SJ_GRID
        global best2, done2
        if past_deadline("stage2")
            @printf("  !! stage-2 deadline (%.1f h) reached after %d/%d points; truncating.
",
                    deadline_h("stage2"), done2, n2)
            break
        end
        p = merge(params, (; J1_meV=j1, sigma_J=sj,
                             gzz=best1.gzz, sigma_gzz=best1.sigma_gzz))
        r = evaluate(p, cuts_disp; label=@sprintf("J1=%.3f sJ=%.2f", j1, sj))
        row!(io2, j1, sj, best1.gzz, best1.sigma_gzz, r.chi2, r.rms, r.scale, r.ok,
             r.n_failed, round(r.seconds; digits=2), r.reg)
        if r.ok && r.chi2 < best2.chi2
            best2 = (; chi2=r.chi2, J1=j1, sigma_J=sj)
        end
        done2 += 1
        @printf("  [%3d/%3d] J1=%.3f sigma_J=%.2f  chi2=%-11.4g %s (%.0f s, %.1f h left)\n",
                done2, n2, j1, sj, r.chi2, r.ok ? " " : "FAIL", r.seconds, remaining_h())
        flush(stdout)
    end
    close(io2)
    @printf("\n  stage 2 best: J1 = %.3f meV, sigma_J = %.2f, chi2_red = %.4g\n",
            best2.J1, best2.sigma_J, best2.chi2)
    println("  -> $p2")
end

# ---------------------------------------------------------------- stage 3

banner("Stage 3 -- consistency: re-fit Gamma at the stage-2 J, coarsely")
println("If (gzz, sigma_gzz) moved, the factorization needs another iteration.\n")

best3 = best1
if past_deadline("stage3")
    println("  !! past the stage-3 deadline; SKIPPING the Gamma recheck.")
else
    io3, p3 = open_csv("stage3_gamma_recheck.csv",
        "gzz,sigma_gzz,J1_meV,sigma_J,chi2_red,rms,scale,ok,seconds")
    gz_c = [best1.gzz + d for d in (-0.2, -0.1, 0.0, 0.1, 0.2)]
    sg_c = [max(0.0, best1.sigma_gzz + d) for d in (-0.3, -0.15, 0.0, 0.15, 0.3)]
    best3 = (; chi2=Inf, gzz=best1.gzz, sigma_gzz=best1.sigma_gzz)
    for gz in gz_c, sg in unique(sg_c)
        global best3
        past_deadline("stage3") && break
        p = merge(params, (; gzz=gz, sigma_gzz=sg, J1_meV=best2.J1, sigma_J=best2.sigma_J))
        r = evaluate(p, cuts_gamma; label=@sprintf("recheck gzz=%.2f", gz))
        row!(io3, gz, sg, best2.J1, best2.sigma_J, r.chi2, r.rms, r.scale, r.ok,
             round(r.seconds; digits=2))
        if r.ok && r.chi2 < best3.chi2
            best3 = (; chi2=r.chi2, gzz=gz, sigma_gzz=sg)
        end
        @printf("  gzz=%.2f sigma_gzz=%.2f  chi2=%-11.4g (%.0f s)\n", gz, sg, r.chi2, r.seconds)
        flush(stdout)
    end
    close(io3)
    @printf("\n  stage 3 best: gzz = %.2f, sigma_gzz = %.2f (stage 1 gave %.2f, %.2f)\n",
            best3.gzz, best3.sigma_gzz, best1.gzz, best1.sigma_gzz)
    moved = abs(best3.gzz - best1.gzz) > 1e-9 || abs(best3.sigma_gzz - best1.sigma_gzz) > 1e-9
    println(moved ?
        "  => the Gamma optimum MOVED under the new J. Iterate: feed these back to stage 2." :
        "  => the Gamma optimum is STABLE. The factorization converged in one pass.")
    println("  -> $p3")
end

# ---------------------------------------------------------------- summary

banner("Stage 4 -- joint local refinement on ALL SIX cuts")
println("""
The sequential factorization optimizes each subset IN TURN, which does not guarantee
that the global chi2 improves. gzz shifts the mode energy at EVERY q, not only at the
zone center, so a gzz chosen from Gamma alone can degrade K and M, and only a joint
step can trade the two off. A smoke test at 12x12x1 showed exactly that failure: Gamma
improved 3.7x while the total chi2 got worse, and the stage-3 "converged" test did not
notice, because it only asked whether Gamma had moved.

So the factorization is used for what it is genuinely good for -- a well-conditioned
starting point that breaks the sigma_J/J1 confound -- and a coordinate-descent sweep on
all six cuts finishes the job.
""")

p_factorized = merge(params, (; gzz=best3.gzz, sigma_gzz=best3.sigma_gzz,
                                J1_meV=best2.J1, sigma_J=best2.sigma_J))
const REFINE_STEPS = (
    (:gzz,        (-0.10, 0.10),  0.0),
    (:sigma_gzz,  (-0.15, 0.15),  0.0),
    (:J1_meV,     (-0.05, 0.05),  0.0),
    (:sigma_J,    (-0.10, 0.10),  0.0),
)
const N_SWEEPS = Int(get(RUN, "refine_sweeps", 2))

p_refined = p_factorized
if past_deadline("stage4")
    println("  !! past the stage-4 deadline; SKIPPING the joint refinement.")
else
    io5, p5 = open_csv("stage4_joint_refinement.csv",
        "sweep,parameter,J1_meV,sigma_J,gzz,sigma_gzz,chi2_red,rms,scale,accepted,seconds")
    r0 = evaluate(p_refined, CUTS_ALL; label="refine-start")
    cur_chi2 = r0.chi2
    row!(io5, 0, "start", p_refined.J1_meV, p_refined.sigma_J, p_refined.gzz,
         p_refined.sigma_gzz, r0.chi2, r0.rms, r0.scale, true, round(r0.seconds; digits=2))
    @printf("  start (factorized): chi2_red = %.5g\n", cur_chi2)
    flush(stdout)
    for sweep in 1:N_SWEEPS
        global p_refined, cur_chi2
        improved = false
        for (key, deltas, floor_val) in REFINE_STEPS
            past_deadline("stage4") && break
            for d in deltas
                past_deadline("stage4") && break
                trial = merge(p_refined,
                    NamedTuple{(key,)}((max(floor_val, getfield(p_refined, key) + d),)))
                getfield(trial, key) == getfield(p_refined, key) && continue
                r = evaluate(trial, CUTS_ALL; label="$key$(d > 0 ? "+" : "")$d")
                acc = r.ok && r.chi2 < cur_chi2
                row!(io5, sweep, String(key), trial.J1_meV, trial.sigma_J, trial.gzz,
                     trial.sigma_gzz, r.chi2, r.rms, r.scale, acc,
                     round(r.seconds; digits=2))
                @printf("  [sweep %d] %-11s %+.3f -> chi2 = %-11.5g %s (%.0f s, %.1f h left)\n",
                        sweep, String(key), d, r.chi2, acc ? "ACCEPT" : "reject",
                        r.seconds, remaining_h())
                flush(stdout)
                if acc
                    p_refined = trial
                    cur_chi2 = r.chi2
                    improved = true
                end
            end
        end
        if !improved
            @printf("  sweep %d found no improvement; the local optimum is reached.\n", sweep)
            break
        end
    end
    close(io5)
    @printf("\n  stage 4: chi2_red %.5g -> %.5g\n", r0.chi2, cur_chi2)
    println("  -> $p5")
end

banner("Summary -- by-eye vs factorized vs jointly refined, on ALL SIX cuts")
p_final = p_refined
io4, p4 = open_csv("summary_all_cuts.csv",
    "variant,J1_meV,sigma_J,gzz,sigma_gzz,chi2_red,rms,scale,field_T,qtag,cut_chi2_red,cut_rms")

totals = Dict{String,Float64}()
for (name, p) in (("by_eye", params), ("factorized", p_factorized), ("refined", p_final))
    remaining_h() < 0.05 && (println("  !! budget exhausted; skipping $name."); break)
    r = evaluate(p, CUTS_ALL; label=name)
    totals[name] = r.chi2
    @printf("  %-12s J1=%.3f sigma_J=%.2f gzz=%.2f sigma_gzz=%.2f  chi2_red = %.5g\n",
            name, p.J1_meV, p.sigma_J, p.gzz, p.sigma_gzz, r.chi2)
    for pc in r.per_cut
        row!(io4, name, p.J1_meV, p.sigma_J, p.gzz, p.sigma_gzz, r.chi2, r.rms, r.scale,
             pc.field_T, pc.qtag, pc.chi2_red, pc.rms)
        @printf("      %-14s %5.1f T  chi2_red = %-11.4g rms = %.4g\n",
                pc.qtag, pc.field_T, pc.chi2_red, pc.rms)
    end
    flush(stdout)
end
close(io4)

# The honest test of the whole exercise. Report a regression loudly rather than
# presenting a "factorized optimum" that is globally worse than where we started.
if haskey(totals, "by_eye")
    be = totals["by_eye"]
    for name in ("factorized", "refined")
        haskey(totals, name) || continue
        d = 100 * (totals[name] - be) / be
        @printf("\n  %s vs by_eye: chi2_red %.5g -> %.5g (%+.1f%%)  %s\n",
                name, be, totals[name], d,
                totals[name] < be ? "IMPROVED" : "*** WORSE THAN THE STARTING POINT ***")
    end
    if haskey(totals, "refined") && totals["refined"] >= be
        println("""
  The refined point is no better than the by-eye parameters on the global chi2. Do NOT
  read that as "the by-eye values are right" -- more likely the four parameters cannot
  simultaneously fit all six cuts, which is a MODEL statement, not a fitting failure.
  The most likely missing ingredient is the one the model does not have: the published
  fits use XXZ with Delta ~ 1.35 while this model is isotropic Heisenberg with all
  anisotropy in the g tensor. Check whether the residual is systematic in q.""")
    end
end

@printf("\nfactorized  (stages 1-3): J1_meV = %.3f, sigma_J = %.2f, gzz = %.2f, sigma_gzz = %.2f\n",
        p_factorized.J1_meV, p_factorized.sigma_J, p_factorized.gzz, p_factorized.sigma_gzz)
@printf("refined     (stage 4)   : J1_meV = %.3f, sigma_J = %.2f, gzz = %.2f, sigma_gzz = %.2f\n",
        p_final.J1_meV, p_final.sigma_J, p_final.gzz, p_final.sigma_gzz)
println("Use the REFINED set. The factorized set is a well-conditioned waypoint, not an answer.")
@printf("total wall clock: %.2f h\n", elapsed_h())
println("tables -> $OUTDIR")
println("\ndone")
