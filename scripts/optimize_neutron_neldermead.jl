#!/usr/bin/env julia
# ONE Nelder-Mead start on the six-cut neutron objective. Run several concurrently, one
# process per start:
#
#   YZGO_START_INDEX=k julia -t 4 --project=. scripts/optimize_neutron_neldermead.jl
#
# Multi-start is naturally parallel at the PROCESS level, which is also the fastest launch
# configuration on this box: intra-process q-threading saturates near 3.4x at 32 threads
# but is ~73% efficient at 4, so 8 processes x 4 threads beats 1 process x 32 threads
# several-fold. It also tests for multiple minima, which a single descent cannot.
#
# This closes the "no optimizer" half of open thread 2. Every evaluation is appended to CSV
# as it completes, so a kill at any point leaves the whole trajectory usable, and a
# wall-clock limit is enforced by Optim so an unattended run cannot overrun.

using Printf, Statistics, LinearAlgebra, Optim, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const LOADED = SV.sv_load_diagnostic_controls(REPO_ROOT,
    "configs/neutron_optimization_controls.toml"; env_var="YZGO_NEUTRON_OPT_CONTROLS")
const CFG = LOADED.diag
const RUN = get(CFG, "run", Dict{String,Any}())
const controls = LOADED.controls

params0, _ = SV.sv_apply_param_overrides(
    SV.sv_load_params(REPO_ROOT, controls).params, RUN)

const CUTS = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
const NREAL = Int(get(RUN, "n_realizations", 4))
const MAXITERS = Int(get(RUN, "minimize_maxiters", 1000))
const RELAX = Int(get(RUN, "relax_attempts", 1))
const TDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "table_subdir",
    "results/feature_tables/sunny_validation/neutron_optimization"))
mkpath(TDIR)

const FREE = Symbol.(String.(get(RUN, "free_parameters",
    ["J1_meV", "sigma_J", "gzz", "sigma_gzz"])))
const LO = Float64.(get(RUN, "lower_bounds", [0.05, 0.0, 2.8, 0.0]))
const HI = Float64.(get(RUN, "upper_bounds", [0.40, 0.60, 4.2, 1.4]))

# sigma_J's upper bound is 0.60 on purpose. It is FRACTIONAL -- J*(1 + sigma_J*randn) --
# so bonds change sign when randn < -1/sigma_J: ~5% of bonds at 0.60 but ~16% at 1.0.
# Beyond ~0.6 this stops being a broadened antiferromagnet and becomes a random-sign bond
# model, which is different physics from the claim being tested.

starts = get(CFG, "starts", Any[])
isempty(starts) && error("No [[starts]] in the config.")
const K = let s = get(ENV, "YZGO_START_INDEX", "0")
    i = something(tryparse(Int, s), 0)
    clamp(i, 0, length(starts) - 1)
end
const START = starts[K + 1]
const TAG = String(get(START, "name", "start$K"))

x0 = [Float64(get(START, String(s), Float64(getproperty(params0, s)))) for s in FREE]
x0 = clamp.(x0, LO, HI)

setp(p, xs) = merge(p, NamedTuple{Tuple(FREE)}(Tuple(xs)))

@printf("start %d (%s), threads %d, %d realizations, %d cuts, %d q/cut\n",
        K, TAG, Threads.nthreads(), NREAL, length(CUTS),
        length(SV.sv_kpm_1d_q_sampler(CUTS[1], controls).qs))
@printf("free: %s\n  lo: %s\n  hi: %s\n  x0: %s\n",
        join(String.(FREE), ", "), string(LO), string(HI), string(round.(x0; digits=4)))
Threads.nthreads() > 8 && @warn "More than 8 threads per process wastes the box in a multi-start run; 4 is the measured sweet spot."

const EVCSV = joinpath(TDIR, "evals_$(TAG).csv")
evio = open(EVCSV, "w")
println(evio, "eval," * join(String.(FREE), ",") * ",chi2_red,rms,scale,ok,n_failed,seconds")
flush(evio)

const T0 = time()
nev = Ref(0)
best = Ref((Inf, copy(x0)))

function objective(x)
    xc = clamp.(x, LO, HI)
    pen = 1e4 * sum(abs, xc .- x)          # soft box; NelderMead is unconstrained
    p = setp(params0, xc)
    t = time()
    local o
    try
        o = SV.sv_neutron_objective(p, controls, CUTS; realizations=0:(NREAL - 1),
            threaded=true, maxiters=MAXITERS, relax_attempts=RELAX, on_failure=:record)
    catch err
        # An unattended optimizer must never die on one bad point. Return a large finite
        # value so the simplex walks away from it instead.
        @printf("    eval %3d FAILED: %s\n", nev[] + 1,
                first(split(sprint(showerror, err), '\n'))[1:min(60, end)])
        flush(stdout)
        nev[] += 1
        println(evio, join(vcat(nev[], round.(xc; digits=5), 1e6, NaN, NaN, false, -1,
                                round(time() - t; digits=1)), ","))
        flush(evio)
        return 1e6 + pen
    end
    nev[] += 1
    f = isfinite(o.chi2_red) ? o.chi2_red : 1e6
    println(evio, join(vcat(nev[], round.(xc; digits=5), round(f; sigdigits=8),
                            round(o.rms; sigdigits=6), round(o.scale; sigdigits=6),
                            o.ok, o.n_failed, round(time() - t; digits=1)), ","))
    flush(evio)
    if f < best[][1]
        best[] = (f, copy(xc))
    end
    @printf("  eval %3d  %s -> chi2 %10.5g  (best %.5g, %.2f h elapsed)\n", nev[],
            join([@sprintf("%s=%.4f", s, v) for (s, v) in zip(FREE, xc)], " "),
            f, best[][1], (time() - T0) / 3600)
    flush(stdout)
    return f + pen
end

const HOURS = Float64(get(RUN, "time_limit_hours", 30.0))
@printf("\nrunning Nelder-Mead, max %d iterations or %.1f h\n\n",
        Int(get(RUN, "max_iterations", 200)), HOURS)

res = Optim.optimize(objective, x0, NelderMead(),
    Optim.Options(iterations = Int(get(RUN, "max_iterations", 200)),
                  time_limit = HOURS * 3600, g_tol = 1e-10, show_trace = false))
close(evio)

xb = clamp.(Optim.minimizer(res), LO, HI)
fb = best[][1]
# Optim's minimizer can differ from the best point actually SEEN when the run stops on a
# limit rather than converging, so report the best seen and trust that.
xseen = best[][2]

open(joinpath(TDIR, "best_$(TAG).csv"), "w") do io
    println(io, "tag,start_index," * join(String.(FREE), ",") *
                ",chi2_red,n_evals,hours,converged,stopped_by")
    println(io, join(vcat(TAG, K, round.(xseen; digits=5), round(fb; sigdigits=8), nev[],
                          round((time() - T0) / 3600; digits=3),
                          Optim.converged(res),
                          Optim.iterations(res) >= Int(get(RUN, "max_iterations", 200)) ?
                              "iterations" : ((time() - T0) >= HOURS * 3600 ? "time" : "converged")), ","))
end

@printf("\n%d evaluations in %.2f h.  %s\n", nev[], (time() - T0) / 3600,
        Optim.converged(res) ? "CONVERGED" : "stopped on a limit")
@printf("best seen chi2_red = %.6g at\n", fb)
for (s, a, b) in zip(FREE, x0, xseen)
    @printf("  %-11s %8.4f -> %8.4f\n", String(s), a, b)
end
maximum(abs.(xb .- xseen)) > 1e-6 &&
    @printf("note: Optim's minimizer differs from the best point seen; reporting best seen.\n")
println("-> $EVCSV")
