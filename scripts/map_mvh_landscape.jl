#!/usr/bin/env julia

# Map the M(H) objective landscape before trying to optimize it.
#
# M(H) is one smooth monotonic digitized curve with no error bars, against a model
# with J1, J2, sigma_J, gzz, sigma_gzz, chi_vv and A_M. That is very likely a
# degenerate fitting problem, and an optimizer would happily return a confident
# point from a flat valley. So: measure which parameters the data can constrain,
# and which combinations they cannot, BEFORE optimizing.
#
# The two parameters the model is linear in are profiled out analytically at every
# grid point by nonnegative least squares (A_M and the Van Vleck slope), so this
# maps the SHAPE residual only. That also makes the whole map insensitive to the
# suspected normalization problem in the digitized data.
#
# Expected degeneracies, stated in advance so they are testable:
#   J1 <-> gzz   through B_sat = S*D_max(J1,J2)/(gzz*mu_B); with A_M profiled out
#                the moment scale carries no information, so J1/gzz should be a
#                near-flat direction.
#   sigma_J <-> sigma_gzz   both round the saturation and add initial slope.
#   J2       likely unconstrained: 0.01 meV entering only via 8(J1+J2) vs 9J1.
#
# Protocol is the one fixed by scripts/check_mvh_convergence.jl: T = 0
# minimize_energy!, field-polarized start, adiabatic continuation, fixed
# realizations (common random numbers) so the objective is deterministic.
#
# Run with THREADS — realizations are independent and this is the whole speedup:
#   julia -t auto --project=. scripts/map_mvh_landscape.jl

using Printf
using Statistics
using LinearAlgebra
using CairoMakie
using Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

const SV = SunnyValidation

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

# Shared script helpers now live in src/sunny_validation.jl (they had been
# copy-pasted across six scripts). These thin aliases keep the local call
# sites unchanged while the logic lives in one place.
const _deepmerge! = SV.sv_deepmerge!

const _repo_path = SV.sv_repo_path

_load_controls(repo_root) = (r = SV.sv_load_diagnostic_controls(repo_root,
        "configs/mvh_landscape_controls.toml"; env_var="SUNNY_MVH_LANDSCAPE_CONTROLS");
    (; r.diag, r.controls, r.diag_path, base=r.base_path))

const _apply_overrides = SV.sv_apply_param_overrides

_sub(d::Dict, k) = get(d, k, Dict{String,Any}())
_set(params, name::AbstractString, v::Real) = merge(params, NamedTuple((Symbol(name) => Float64(v),)))

"Scan range for a parameter: explicit absolute range if given, else +-span around the centre."
function _range_for(name::AbstractString, centre::Real, run::Dict, sub::Dict)
    abs_ranges = _sub(sub, "absolute_ranges")
    if haskey(abs_ranges, name)
        r = abs_ranges[name]
        return (Float64(r[1]), Float64(r[2]))
    end
    span = Float64(get(sub, "relative_span", get(run, "relative_span", 0.6)))
    c = Float64(centre)
    return (max(0.0, c * (1 - span)), c * (1 + span))
end

# -----------------------------------------------------------------------------
# Objective wrapper
# -----------------------------------------------------------------------------

struct Obj
    controls::Dict
    Bs::Vector{Float64}
    M_exp::Vector{Float64}
    cell_size
    seed_dims
    realizations
    maxiters::Int
    threaded::Bool
end

function (o::Obj)(params)
    return SV.sv_mvh_objective(params, o.controls, o.Bs, o.M_exp;
        cell_size=o.cell_size, seed_dims=o.seed_dims, realizations=o.realizations,
        maxiters=o.maxiters, threaded=o.threaded)
end

_write_csv(path, header, rows) = SV.sv_write_rows_csv(path, rows; header)

_hdr(rows) = String.(collect(keys(rows[1])))

# -----------------------------------------------------------------------------
# Reproducibility floor: how much does the objective move if only the disorder
# realization SET changes? Differences smaller than this are meaningless, so this
# is the yardstick for calling a direction flat.
# -----------------------------------------------------------------------------

function reproducibility_floor(o::Obj, params)
    nR = length(o.realizations)
    alt = (last(o.realizations) + 1):(last(o.realizations) + nR)
    o2 = Obj(o.controls, o.Bs, o.M_exp, o.cell_size, o.seed_dims, alt, o.maxiters, o.threaded)
    a, b = o(params), o2(params)
    floor = abs(a.rms - b.rms)
    println("\n================ reproducibility floor ================")
    @printf("  realizations %s : rms = %.5f  A_M = %.4f  chi_vv = %.4f  (%.1f s)\n",
            string(o.realizations), a.rms, a.A_M, a.chi_vv, a.seconds)
    @printf("  realizations %s : rms = %.5f  A_M = %.4f  chi_vv = %.4f  (%.1f s)\n",
            string(alt), b.rms, b.A_M, b.chi_vv, b.seconds)
    @printf("  |d rms| from changing the realization set = %.5f uB\n", floor)
    println("  Any variation in the maps below this is not meaningful.")
    return (; floor, a, b, alt)
end

# -----------------------------------------------------------------------------
# 1D scans
# -----------------------------------------------------------------------------

function scan_1d(o::Obj, params0, run::Dict, out, floor::Float64)
    s = _sub(run, "scan1d")
    Bool(get(s, "enabled", true)) || return nothing
    names = String.(get(s, "parameters", ["J1_meV", "sigma_J", "gzz", "sigma_gzz", "J2_meV"]))
    n = Int(get(s, "n_points", 11))

    println("\n================ 1D sensitivity scans ================")
    rows = NamedTuple[]
    summary = NamedTuple[]
    for name in names
        c = Float64(getproperty(params0, Symbol(name)))
        lo, hi = _range_for(name, c, run, s)
        xs = collect(range(lo, hi; length=n))
        @printf("\n  %s : %.4g -> %.4g, centre %.4g\n", name, lo, hi, c)
        @printf("  %12s %10s %10s %10s\n", name, "rms", "A_M", "chi_vv")
        best = (rms=Inf, x=NaN)
        for x in xs
            r = o(_set(params0, name, x))
            push!(rows, (; parameter=name, value=x, rms=r.rms, max_abs=r.max_abs,
                A_M=r.A_M, chi_vv=r.chi_vv, nan_fraction=r.nan_fraction))
            r.rms < best.rms && (best = (rms=r.rms, x=x))
            @printf("  %12.5g %10.5f %10.4f %10.4f\n", x, r.rms, r.A_M, r.chi_vv)
        end
        sel = [x for x in rows if x.parameter == name]
        rmin, rmax = minimum(x.rms for x in sel), maximum(x.rms for x in sel)
        # Flat interval: the span of parameter values whose rms is within the
        # reproducibility floor of the best point. That is the honest uncertainty.
        within = [x.value for x in sel if x.rms <= rmin + floor]
        w = isempty(within) ? NaN : (maximum(within) - minimum(within))
        constrained = (rmax - rmin) > 3 * floor
        push!(summary, (; parameter=name, centre=c, best=best.x, rms_min=rmin,
            rms_range=rmax - rmin, flat_span=w,
            flat_span_frac_of_range=w / (hi - lo), constrained))
        @printf("  -> best %.5g, rms range %.5f, within-floor span %.4g (%.0f%% of scan)  %s\n",
                best.x, rmax - rmin, w, 100 * w / (hi - lo),
                constrained ? "CONSTRAINED" : "*** FLAT / UNCONSTRAINED ***")
    end

    println("\n  summary:")
    @printf("  %-12s %10s %10s %12s %10s  %s\n",
            "parameter", "centre", "best", "rms range", "flat span", "verdict")
    for x in summary
        @printf("  %-12s %10.4g %10.4g %12.5f %10.4g  %s\n", x.parameter, x.centre,
                x.best, x.rms_range, x.flat_span,
                x.constrained ? "constrained" : "FLAT")
    end

    _write_csv(joinpath(out, "mvh_landscape_scan1d.csv"), _hdr(rows), rows)
    _write_csv(joinpath(out, "mvh_landscape_scan1d_summary.csv"), _hdr(summary), summary)
    return (; rows, summary, names)
end

# -----------------------------------------------------------------------------
# 2D maps over the suspected degenerate pairs
# -----------------------------------------------------------------------------

function scan_2d(o::Obj, params0, run::Dict, out, floor::Float64)
    s = _sub(run, "scan2d")
    Bool(get(s, "enabled", true)) || return nothing
    pairs = get(s, "pairs", [["sigma_J", "sigma_gzz"], ["J1_meV", "gzz"]])
    n = Int(get(s, "n_points", 11))

    println("\n================ 2D degeneracy maps ================")
    maps = []
    rows = NamedTuple[]
    for pr in pairs
        n1, n2 = String(pr[1]), String(pr[2])
        c1 = Float64(getproperty(params0, Symbol(n1)))
        c2 = Float64(getproperty(params0, Symbol(n2)))
        l1, h1 = _range_for(n1, c1, run, s)
        l2, h2 = _range_for(n2, c2, run, s)
        xs, ys = collect(range(l1, h1; length=n)), collect(range(l2, h2; length=n))
        Z = fill(NaN, n, n)
        @printf("\n  %s x %s : [%.4g, %.4g] x [%.4g, %.4g], %dx%d\n",
                n1, n2, l1, h1, l2, h2, n, n)
        t = @elapsed for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
            p = _set(_set(params0, n1, x), n2, y)
            r = o(p)
            Z[i, j] = r.rms
            push!(rows, (; par1=n1, par2=n2, val1=x, val2=y, rms=r.rms,
                A_M=r.A_M, chi_vv=r.chi_vv))
        end
        idx = argmin(replace(Z, NaN => Inf))
        zmin = Z[idx]
        nflat = count(z -> isfinite(z) && z <= zmin + floor, Z)
        @printf("  %.1f s   best (%s=%.4g, %s=%.4g) rms=%.5f\n", t,
                n1, xs[idx[1]], n2, ys[idx[2]], zmin)
        @printf("  %d of %d grid cells lie within the reproducibility floor of the best\n",
                nflat, n * n)
        @printf("  -> %s\n", nflat > 0.2 * n * n ?
                "*** BROAD DEGENERATE VALLEY — this pair is not separately determined ***" :
                "localized minimum, pair is separable")
        push!(maps, (; n1, n2, xs, ys, Z, zmin, best=(xs[idx[1]], ys[idx[2]]), nflat))
    end
    _write_csv(joinpath(out, "mvh_landscape_scan2d.csv"), _hdr(rows), rows)
    return (; maps, rows)
end

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------

function make_figure(path, s1, s2, floor, params0)
    mkpath(dirname(path))
    nmaps = s2 === nothing ? 0 : length(s2.maps)
    fig = Figure(size=(560 * max(2, 1 + nmaps), 900))
    pal = [:royalblue, :darkorange, :seagreen, :orchid, :grey40, :firebrick]

    # A: all 1D scans on a normalized axis, so steep and flat are directly comparable.
    ax = Axis(fig[1, 1]; xlabel="parameter, normalized over its scan range",
        ylabel="rms residual (μB/Yb)",
        title="A. 1D sensitivity — flat curves are unconstrained directions")
    if s1 !== nothing
        for (k, name) in enumerate(s1.names)
            sel = [x for x in s1.rows if x.parameter == name]
            isempty(sel) && continue
            vs = [x.value for x in sel]
            lo, hi = minimum(vs), maximum(vs)
            xn = hi > lo ? (vs .- lo) ./ (hi - lo) : fill(0.5, length(vs))
            scatterlines!(ax, xn, [x.rms for x in sel]; color=pal[mod1(k, 6)],
                label=@sprintf("%s [%.3g, %.3g]", name, lo, hi))
        end
        rmin = minimum(x.rms for x in s1.rows)
        hlines!(ax, [rmin + floor]; color=:black, linestyle=:dash,
            label=@sprintf("best + reproducibility floor (%.4f)", floor))
        axislegend(ax; position=:rt, labelsize=8, framevisible=false)
    end

    # B: how the profiled-out linear parameters move along the scans. If A_M has
    # to swing wildly to keep the fit, that IS the degeneracy.
    ax = Axis(fig[2, 1]; xlabel="parameter, normalized over its scan range",
        ylabel="fitted A_M", title="B. Profiled-out amplitude along each scan")
    if s1 !== nothing
        for (k, name) in enumerate(s1.names)
            sel = [x for x in s1.rows if x.parameter == name]
            isempty(sel) && continue
            vs = [x.value for x in sel]
            lo, hi = minimum(vs), maximum(vs)
            xn = hi > lo ? (vs .- lo) ./ (hi - lo) : fill(0.5, length(vs))
            scatterlines!(ax, xn, [x.A_M for x in sel]; color=pal[mod1(k, 6)], label=name)
        end
        axislegend(ax; position=:rt, labelsize=8, framevisible=false)
    end

    # C, D: the 2D maps. The floor contour is the confidence region.
    for (m, mp) in enumerate(s2 === nothing ? [] : s2.maps)
        for (row, what) in ((1, :rms), (2, :zoom))
            ax = Axis(fig[row, 1 + m]; xlabel=mp.n1, ylabel=mp.n2,
                title=row == 1 ?
                    @sprintf("%s. rms over (%s, %s)", "CD"[min(m, 2)], mp.n1, mp.n2) :
                    "within reproducibility floor of best (white)")
            if what == :rms
                hm = heatmap!(ax, mp.xs, mp.ys, mp.Z; colormap=:viridis)
                contour!(ax, mp.xs, mp.ys, mp.Z;
                    levels=[mp.zmin + floor, mp.zmin + 2 * floor, mp.zmin + 4 * floor],
                    color=:white, linewidth=1.2)
                Colorbar(fig[row, 1 + m, Right()], hm; width=12, ticklabelsize=8,
                    tickformat=vs -> [@sprintf("%.3f", v) for v in vs])
            else
                mask = [isfinite(z) && z <= mp.zmin + floor ? 1.0 : 0.0 for z in mp.Z]
                heatmap!(ax, mp.xs, mp.ys, mask; colormap=:grays)
            end
            scatter!(ax, [mp.best[1]], [mp.best[2]]; color=:red, marker=:xcross, markersize=14)
            scatter!(ax, [Float64(getproperty(params0, Symbol(mp.n1)))],
                        [Float64(getproperty(params0, Symbol(mp.n2)))];
                color=:cyan, marker=:circle, markersize=10)
        end
    end

    Label(fig[0, :], "M(H) objective landscape — shape only, A_M and χ_vv profiled out (red × = best, cyan ○ = starting point)";
        fontsize=15, font=:bold)
    save(path, fig)
    return path
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

function main()
    (; diag, controls, diag_path, base) = _load_controls(REPO_ROOT)
    run = get(diag, "run", Dict{String,Any}())
    (; params, path) = SV.sv_load_params(REPO_ROOT, controls)
    params0, overrides = _apply_overrides(params, run)

    @info "M(H) objective landscape" sunny_version=SV.sv_try_pkgversion(Sunny) params_path=path diag_path base
    if !isempty(overrides)
        println("\nCentre point [run.param_overrides]:")
        for s in overrides
            println("  ", s)
        end
    end

    nB = Int(get(run, "nB", 18))
    Bs = collect(range(Float64(get(run, "Bmin_T", 0.2)), Float64(get(run, "Bmax_T", 6.8)); length=nB))
    M_exp = SV.sv_mvh_target(REPO_ROOT, controls, Bs)
    nR = Int(get(run, "n_realizations", 16))
    off = Int(get(run, "realization_offset", 0))

    if Threads.nthreads() == 1
        @warn "Running single-threaded. Realizations are independent — relaunch with `julia -t auto` for a large speedup."
    end
    @printf("\nThreads: %d.  Realizations: %d (offset %d).  Fields: %d over %.3g-%.3g T.\n",
            Threads.nthreads(), nR, off, nB, first(Bs), last(Bs))
    @printf("Cell %s, seed %s, maxiters %d.  Experiment points valid: %d/%d.\n",
            SV.sv_cell_label(get(run, "cell_size", [12, 12, 1])),
            SV.sv_cell_label(get(run, "seed_dims", [3, 3, 1])),
            Int(get(run, "minimize_maxiters", 30_000)), count(isfinite, M_exp), nB)

    o = Obj(controls, Bs, M_exp,
            get(run, "cell_size", [12, 12, 1]), get(run, "seed_dims", [3, 3, 1]),
            off:(off + nR - 1), Int(get(run, "minimize_maxiters", 30_000)),
            Bool(get(run, "threaded", true)))

    out_tab = SV.sv_repo_path(REPO_ROOT, controls["paths"]["table_subdir"])
    out_fig = SV.sv_repo_path(REPO_ROOT, controls["paths"]["figure_subdir"])
    mkpath(out_tab); mkpath(out_fig)

    rf = reproducibility_floor(o, params0)
    total = @elapsed begin
        s1 = scan_1d(o, params0, run, out_tab, rf.floor)
        s2 = scan_2d(o, params0, run, out_tab, rf.floor)
    end
    fig_path = make_figure(joinpath(out_fig, "mvh_landscape.png"), s1, s2, rf.floor, params0)

    println("\n================ what M(H) can and cannot constrain ================")
    if s1 !== nothing
        con = [x.parameter for x in s1.summary if x.constrained]
        flat = [x.parameter for x in s1.summary if !x.constrained]
        println("  constrained by M(H) shape : ", isempty(con) ? "(none)" : join(con, ", "))
        println("  flat / unconstrained      : ", isempty(flat) ? "(none)" : join(flat, ", "))
        isempty(flat) || println("  Flat parameters must be pinned from another observable, not fitted here.")
    end
    if s2 !== nothing
        for mp in s2.maps
            frac = mp.nflat / length(mp.Z)
            @printf("  (%s, %s): %.0f%% of the grid is within the floor -> %s\n",
                    mp.n1, mp.n2, 100 * frac,
                    frac > 0.2 ? "degenerate pair" : "separable")
        end
    end
    @printf("\nScan compute: %.1f s\n", total)
    println("Tables: ", out_tab)
    println("Figure: ", fig_path)
    return (; s1, s2, floor=rf.floor, fig_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
