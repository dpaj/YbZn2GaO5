#!/usr/bin/env julia

# Convergence and protocol diagnostics for the large-cell M(H) model.
#
# Purpose: decide, on measured numbers rather than intuition, what cell size,
# how many disorder realizations, and which relaxation protocol give a converged
# and trustworthy M(H). That is the prerequisite for using M(H) as an optimizer
# objective, and later for co-optimizing it against the neutron spectra.
#
# Everything runs at T = 0 with minimize_energy!, because the M(H,T) study showed
# classical finite temperature changes M(H) by less than the model-experiment
# residual at 0.42 K (RMS 0.0182 vs 0.0176 uB). Dropping the thermostat makes the
# objective deterministic given a realization, which is what an optimizer needs.
#
# Five studies, each independently toggleable:
#
#   1. scaling         How realization scatter falls with cell size, and whether
#                      the distribution over realizations is Gaussian or
#                      heavy-tailed. A heavy tail means rare weakly-coupled
#                      regions dominate, mean +- std understates the error, and
#                      many more realizations are needed.
#
#   2. initialization  field_polarized vs random vs annealed, at fixed
#                      realizations. If M(H) depends on the starting
#                      configuration then minimize_energy! is trapping in
#                      metastable states and the "ground state" is protocol
#                      dependent — which no cell size fixes.
#
#   3. hysteresis      Sweep the field up, then back down. Any difference is
#                      direct evidence of metastability and domain physics.
#
#   4. structure       Transverse static structure factor and energy per site of
#                      the relaxed state versus cell size. This is the direct
#                      answer to "is the box big enough to hold the texture":
#                      if the structure-factor width is L independent, the
#                      correlation length is resolved.
#
#   5. commensurability  The 120-degree order is a three-sublattice structure, so
#                      cells divisible by 3 accommodate it and others frustrate
#                      it. 16x16 is a deliberate control against 12x12 and 24x24.
#
# Run with:
#   julia --project=. scripts/check_mvh_convergence.jl
#   julia --project=. scripts/check_mvh_convergence.jl configs/my_variant.toml

using Printf
using Statistics
using LinearAlgebra
using Random
using CairoMakie
using Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

const SV = SunnyValidation
const KB_MEV_PER_K = 0.08617333262

# -----------------------------------------------------------------------------
# Config plumbing
# -----------------------------------------------------------------------------

# Shared script helpers now live in src/sunny_validation.jl (they had been
# copy-pasted across six scripts). These thin aliases keep the local call
# sites unchanged while the logic lives in one place.
const _deepmerge! = SV.sv_deepmerge!

const _repo_path = SV.sv_repo_path

_load_controls(repo_root) = (r = SV.sv_load_diagnostic_controls(repo_root,
        "configs/mvh_convergence_controls.toml"; env_var="SUNNY_MVH_CONV_CONTROLS");
    (; r.diag, r.controls, r.diag_path, base=r.base_path))

const _apply_param_overrides = SV.sv_apply_param_overrides

_sub(run::Dict, name) = get(run, name, Dict{String,Any}())
_on(run::Dict, name) = Bool(get(_sub(run, name), "enabled", true))

# -----------------------------------------------------------------------------
# Core: one T = 0 field sweep
# -----------------------------------------------------------------------------

struct Ctx
    params
    controls::Dict
    run::Dict
    uhat
    units
    moment_sign::Float64
    spin_S::Float64
    seed_dims
    maxiters::Int
    kT::Float64
end

function _ctx(params, controls, run)
    Ctx(params, controls, run,
        SV.sv_field_direction(controls), SV.sv_units(),
        Float64(get(get(controls, "largecell", Dict{String,Any}()), "moment_sign", -1.0)),
        Float64(controls["common"]["spin_S"]),
        get(run, "seed_dims", [3, 3, 1]),
        Int(get(run, "minimize_maxiters", 30_000)),
        KB_MEV_PER_K * Float64(get(run, "temperature_K", 0.42)))
end

function _fresh(c::Ctx, cell, realization; field_T=0.0)
    lc = get(c.controls, "largecell", Dict{String,Any}())
    return SV.sv_build_supercell_system(c.params, c.controls;
        component=:dispersive, cell_size=cell, seed_dims=c.seed_dims,
        field_T=field_T, realization=realization,
        include_exchange=Bool(get(lc, "include_exchange_disorder", true)),
        include_gzz=Bool(get(lc, "include_gzz_disorder", true))).sys
end

function _initialize!(c::Ctx, sys, init::AbstractString, realization::Int; field_T=1.0)
    if init == "field_polarized"
        SV.sv_polarize_along_field!(sys, c.uhat; field_T)
    elseif init == "random"
        Random.seed!(Int(c.controls["common"]["seed"]) + 101 + 7919 * realization)
        randomize_spins!(sys)
    elseif init == "annealed"
        an = _sub(c.run, "annealing")
        Random.seed!(Int(c.controls["common"]["seed"]) + 202 + 7919 * realization)
        randomize_spins!(sys)
        SV.sv_anneal_to_ground_state!(sys;
            kT_high=Float64(get(an, "kT_high_meV", 0.3)),
            kT_low=Float64(get(an, "kT_low_meV", c.kT)),
            n_stages=Int(get(an, "n_stages", 8)),
            n_steps=Int(get(an, "n_steps", 400)),
            dt=Float64(get(an, "dt", 0.02)),
            damping=Float64(get(an, "damping", 0.1)),
            maxiters=c.maxiters)
    else
        error("Unknown initialization '$init'. Use field_polarized, random, or annealed.")
    end
    return sys
end

"""
T = 0 field sweep with adiabatic continuation. `Bs` may be non-monotonic, which is
how the up-then-down hysteresis loop is driven.
"""
function _sweep(c::Ctx, cell, realization::Int, init::AbstractString, Bs)
    sys = _fresh(c, cell, realization; field_T=0.0)
    _initialize!(c, sys, init, realization)
    M = zeros(Float64, length(Bs))
    E = zeros(Float64, length(Bs))
    for (i, B) in enumerate(Bs)
        SV.sv_set_field_T!(sys, c.uhat, c.units, B)
        minimize_energy!(sys; maxiters=c.maxiters)
        M[i] = c.moment_sign * SV.sv_m_parallel_uB_per_site(sys, c.uhat)
        E[i] = energy_per_site(sys)
    end
    return (; M, E, M_sat=SV.sv_m_sat_uB_per_site(sys, c.uhat, c.spin_S), sys)
end

# -----------------------------------------------------------------------------
# Distribution summary. Skew and excess kurtosis are the tail diagnostics: a
# Gaussian has 0 for both, and a rare-region-dominated observable is strongly
# positively skewed with heavy tails.
# -----------------------------------------------------------------------------

function _dist(x::AbstractVector{<:Real})
    v = collect(skipmissing(x))
    v = v[isfinite.(v)]
    n = length(v)
    n == 0 && return (; n=0, mean=NaN, std=NaN, sem=NaN, skew=NaN, exkurt=NaN,
                        min=NaN, max=NaN, median=NaN, p90_over_sigma=NaN)
    mu = mean(v)
    sd = n > 1 ? std(v) : 0.0
    z = sd > 0 ? (v .- mu) ./ sd : zeros(n)
    return (; n, mean=mu, std=sd, sem=n > 1 ? sd / sqrt(n) : 0.0,
              skew=n > 2 ? mean(z .^ 3) : NaN,
              exkurt=n > 3 ? mean(z .^ 4) - 3 : NaN,
              min=minimum(v), max=maximum(v), median=median(v),
              p90_over_sigma=sd > 0 ? (quantile(v, 0.9) - mu) / sd : NaN)
end

# -----------------------------------------------------------------------------
# CSV helper
# -----------------------------------------------------------------------------

_write_csv(path, header, rows) = SV.sv_write_rows_csv(path, rows; header)

_hdr(rows) = String.(collect(keys(rows[1])))

# -----------------------------------------------------------------------------
# Study 1: realization scatter versus cell size, and the tail shape
# -----------------------------------------------------------------------------

function study_scaling(c::Ctx, out)
    s = _sub(c.run, "scaling")
    cells = get(s, "cells", [[12, 12, 1], [24, 24, 1], [36, 36, 1]])
    nr = Int(get(s, "n_realizations", 16))
    Bs = Float64.(get(s, "fields_T", [0.5, 1.0, 2.0, 3.0, 5.0, 7.0]))
    init = String(get(s, "initialization", "field_polarized"))

    println("\n================ 1. realization scatter vs cell size ================")
    @printf("%d realizations, fields %s T, init=%s\n", nr,
            join(string.(Bs), "/"), init)

    curves = NamedTuple[]
    for cell in cells
        label = SV.sv_cell_label(cell)
        nsites = prod(SV.sv_tuple3(cell))
        t = @elapsed for r in 0:(nr - 1)
            sw = _sweep(c, cell, r, init, Bs)
            for (i, B) in enumerate(Bs)
                push!(curves, (; cell=label, nsites, realization=r, field_T=B,
                    M=sw.M[i], E_per_site=sw.E[i], M_sat=sw.M_sat))
            end
        end
        @printf("  %-9s %5d sites  %6.1f s\n", label, nsites, t)
    end

    stats = NamedTuple[]
    for cell in cells, B in Bs
        label = SV.sv_cell_label(cell)
        vals = [x.M for x in curves if x.cell == label && x.field_T == B]
        d = _dist(vals)
        push!(stats, (; cell=label, nsites=prod(SV.sv_tuple3(cell)), field_T=B,
            n=d.n, mean=d.mean, std=d.std, sem=d.sem, median=d.median,
            skew=d.skew, excess_kurtosis=d.exkurt, min=d.min, max=d.max,
            p90_minus_mean_over_sigma=d.p90_over_sigma))
    end

    println("\n  scatter (std over realizations) of M, uB/site:")
    @printf("  %8s", "B(T)")
    for cell in cells
        @printf("%12s", SV.sv_cell_label(cell))
    end
    println("   ratio(first/last)   1/sqrt(N) expects")
    for B in Bs
        @printf("  %8.2f", B)
        sds = Float64[]
        for cell in cells
            st = only([x for x in stats if x.cell == SV.sv_cell_label(cell) && x.field_T == B])
            push!(sds, st.std)
            @printf("%12.5f", st.std)
        end
        exp_ratio = sqrt(prod(SV.sv_tuple3(cells[end])) / prod(SV.sv_tuple3(cells[1])))
        @printf("%14.2f%16.2f\n", sds[end] > 0 ? sds[1] / sds[end] : NaN, exp_ratio)
    end

    println("\n  tail shape (Gaussian => skew 0, excess kurtosis 0):")
    @printf("  %-9s %8s %8s %8s %10s   %s\n", "cell", "B(T)", "skew", "exkurt", "max-mean", "verdict")
    for st in stats
        heavy = (isfinite(st.skew) && abs(st.skew) > 1.0) ||
                (isfinite(st.excess_kurtosis) && st.excess_kurtosis > 1.5)
        @printf("  %-9s %8.2f %8.2f %8.2f %10.4f   %s\n", st.cell, st.field_T,
                st.skew, st.excess_kurtosis, st.max - st.mean,
                heavy ? "HEAVY TAIL — needs more realizations" : "consistent with Gaussian")
    end

    _write_csv(joinpath(out, "mvh_convergence_scaling.csv"), _hdr(curves), curves)
    _write_csv(joinpath(out, "mvh_convergence_scaling_stats.csv"), _hdr(stats), stats)
    return (; curves, stats, cells, Bs)
end

# -----------------------------------------------------------------------------
# Study 2: does the answer depend on how we start?
# -----------------------------------------------------------------------------

function study_initialization(c::Ctx, out)
    s = _sub(c.run, "initialization")
    cells = get(s, "cells", [[12, 12, 1], [36, 36, 1]])
    inits = String.(get(s, "inits", ["field_polarized", "random", "annealed"]))
    nr = Int(get(s, "n_realizations", 4))
    Bs = Float64.(get(s, "fields_T", [0.5, 1.0, 2.0, 3.0, 5.0, 7.0]))

    println("\n================ 2. initialization sensitivity ================")
    println("  If M(H) depends on the start, minimize_energy! is trapping and the")
    println("  ground state is protocol dependent — no cell size fixes that.")

    rows = NamedTuple[]
    for cell in cells, init in inits
        label = SV.sv_cell_label(cell)
        t = @elapsed for r in 0:(nr - 1)
            sw = _sweep(c, cell, r, init, Bs)
            for (i, B) in enumerate(Bs)
                push!(rows, (; cell=label, init, realization=r, field_T=B,
                    M=sw.M[i], E_per_site=sw.E[i]))
            end
        end
        @printf("  %-9s %-16s %6.1f s\n", label, init, t)
    end

    println("\n  mean M and mean E/site by initialization (mean over realizations):")
    summ = NamedTuple[]
    for cell in cells
        label = SV.sv_cell_label(cell)
        println("\n  cell = $label")
        @printf("  %8s", "B(T)")
        for init in inits
            @printf("%14s", init)
        end
        @printf("%12s%14s\n", "spread(M)", "best E/site")
        for B in Bs
            @printf("  %8.2f", B)
            ms = Float64[]
            es = Float64[]
            for init in inits
                v = [x.M for x in rows if x.cell == label && x.init == init && x.field_T == B]
                e = [x.E_per_site for x in rows if x.cell == label && x.init == init && x.field_T == B]
                push!(ms, mean(v)); push!(es, mean(e))
                @printf("%14.5f", mean(v))
            end
            spread = maximum(ms) - minimum(ms)
            @printf("%12.5f%14s\n", spread, inits[argmin(es)])
            push!(summ, (; cell=label, field_T=B, spread_M=spread,
                lowest_energy_init=inits[argmin(es)],
                E_spread_meV=maximum(es) - minimum(es)))
        end
    end

    worst = maximum(x.spread_M for x in summ)
    println()
    if worst < 0.005
        @printf("  VERDICT: initialization independent to %.4f uB — protocol is safe.\n", worst)
    else
        @printf("  VERDICT: initialization matters, up to %.4f uB. minimize_energy! is\n", worst)
        println("  finding different local minima. Prefer the lowest-energy protocol and")
        println("  treat single-shot field_polarized results with caution.")
    end

    _write_csv(joinpath(out, "mvh_convergence_initialization.csv"), _hdr(rows), rows)
    _write_csv(joinpath(out, "mvh_convergence_initialization_summary.csv"), _hdr(summ), summ)
    return (; rows, summ, cells, inits, Bs, worst)
end

# -----------------------------------------------------------------------------
# Study 3: field-history dependence
# -----------------------------------------------------------------------------

function study_hysteresis(c::Ctx, out)
    s = _sub(c.run, "hysteresis")
    cells = get(s, "cells", [[12, 12, 1], [36, 36, 1]])
    nr = Int(get(s, "n_realizations", 4))
    nB = Int(get(s, "nB", 15))
    Bmax = Float64(get(s, "Bmax_T", 7.0))
    init = String(get(s, "initialization", "field_polarized"))

    up = collect(range(0.0, Bmax; length=nB))
    loop = vcat(up, reverse(up)[2:end])

    println("\n================ 3. field-history (hysteresis) ================")
    @printf("  up then down, %d points each way, %d realizations, init=%s\n", nB, nr, init)

    rows = NamedTuple[]
    for cell in cells
        label = SV.sv_cell_label(cell)
        t = @elapsed for r in 0:(nr - 1)
            sw = _sweep(c, cell, r, init, loop)
            for (i, B) in enumerate(loop)
                push!(rows, (; cell=label, realization=r, step=i, field_T=B,
                    branch=i <= nB ? "up" : "down", M=sw.M[i], E_per_site=sw.E[i]))
            end
        end
        @printf("  %-9s %6.1f s\n", label, t)
    end

    println("\n  |M_up - M_down| at matched field (mean over realizations):")
    println("  B = 0 is reported separately: there M vanishes by symmetry, so the gap only")
    println("  records which frozen texture was landed in and carries no weight in a fit.")
    println("  The experimental curve also starts at 0.022 T.")
    summ = NamedTuple[]
    worst_nonzero = Dict{String,Float64}()
    for cell in cells
        label = SV.sv_cell_label(cell)
        worst_all, worst_pos = 0.0, 0.0
        for (i, B) in enumerate(up)
            j = length(loop) - i + 1   # matching point on the down branch
            u = [x.M for x in rows if x.cell == label && x.step == i]
            d = [x.M for x in rows if x.cell == label && x.step == j]
            gap = abs(mean(u) - mean(d))
            worst_all = max(worst_all, gap)
            B > 1e-9 && (worst_pos = max(worst_pos, gap))
            push!(summ, (; cell=label, field_T=B, M_up=mean(u), M_down=mean(d), gap))
        end
        worst_nonzero[label] = worst_pos
        @printf("  %-9s worst gap: %.5f at B=0, %.5f for B>0  %s\n", label,
                worst_all, worst_pos,
                worst_pos < 0.005 ? "(reversible where it matters)" :
                                    "(HYSTERETIC — metastability present)")
    end

    _write_csv(joinpath(out, "mvh_convergence_hysteresis.csv"), _hdr(rows), rows)
    _write_csv(joinpath(out, "mvh_convergence_hysteresis_summary.csv"), _hdr(summ), summ)
    return (; rows, summ, cells, up, worst_nonzero)
end

# -----------------------------------------------------------------------------
# Study 4: is the box big enough to hold the texture?
# -----------------------------------------------------------------------------

function study_structure(c::Ctx, out)
    s = _sub(c.run, "structure")
    cells = get(s, "cells", [[12, 12, 1], [24, 24, 1], [36, 36, 1]])
    nr = Int(get(s, "n_realizations", 4))
    Bs = Float64.(get(s, "fields_T", [0.0, 1.0, 3.0, 5.0]))
    init = String(get(s, "initialization", "field_polarized"))

    println("\n================ 4. ground-state texture vs cell size ================")
    println("  width_rlu independent of L => correlation length resolved, small cell OK.")
    println("  width_rlu tracking 1/L     => resolution limited, xi exceeds the cell.")

    rows = NamedTuple[]
    for cell in cells
        label = SV.sv_cell_label(cell)
        L = SV.sv_tuple3(cell)[1]
        t = @elapsed for r in 0:(nr - 1)
            for B in Bs
                sys = _fresh(c, cell, r; field_T=B)
                _initialize!(c, sys, init, r; field_T=max(B, 1.0))
                SV.sv_set_field_T!(sys, c.uhat, c.units, B)
                minimize_energy!(sys; maxiters=c.maxiters)
                sf = SV.sv_transverse_structure_factor(sys, c.uhat)
                push!(rows, (; cell=label, L, realization=r, field_T=B,
                    E_per_site=energy_per_site(sys),
                    M=c.moment_sign * SV.sv_m_parallel_uB_per_site(sys, c.uhat),
                    peak_q1=sf.peak_q[1], peak_q2=sf.peak_q[2],
                    peak_at_K=(abs(sf.peak_q[1] - 1 / 3) < 0.5 / L &&
                               abs(sf.peak_q[2] - 1 / 3) < 0.5 / L),
                    K_fraction=sf.K_fraction,
                    peak_fraction=sf.peak_fraction, width_rlu=sf.width_rlu,
                    xi_estimate_cells=sf.xi_estimate_cells,
                    inv_L=1 / L, width_times_L=sf.width_rlu * L))
            end
        end
        @printf("  %-9s %6.1f s\n", label, t)
    end

    println("\n  transverse structure factor and energy per site:")
    println("  S_K = weight at the 120-degree ordering vector K=(1/3,1/3); peak@K = is the")
    println("  global maximum there. Low S_K with a broad width means no coherent 120-degree")
    println("  texture at all, which is a physics result, not a finite-size artifact.")
    @printf("  %8s %-9s %14s %7s %8s %10s %9s %9s %12s\n",
            "B(T)", "cell", "peak_q", "peak@K", "S_K", "peak_frac", "width", "width*L", "E/site(meV)")
    for B in Bs
        for cell in cells
            label = SV.sv_cell_label(cell)
            sel = [x for x in rows if x.cell == label && x.field_T == B]
            isempty(sel) && continue
            nK = count(x -> x.peak_at_K, sel)
            @printf("  %8.2f %-9s %6.3f,%-7.3f %3d/%-3d %8.4f %10.4f %9.5f %9.3f %12.6f\n",
                    B, label, mean(x.peak_q1 for x in sel), mean(x.peak_q2 for x in sel),
                    nK, length(sel), mean(x.K_fraction for x in sel),
                    mean(x.peak_fraction for x in sel),
                    mean(x.width_rlu for x in sel),
                    mean(x.width_times_L for x in sel), mean(x.E_per_site for x in sel))
        end
        println()
    end

    println("  Reading: if width*L is roughly CONSTANT across cells the peak is")
    println("  resolution limited (xi >= L, cell too small). If width_rlu is roughly")
    println("  constant instead, xi is resolved and the small cell is adequate.")

    _write_csv(joinpath(out, "mvh_convergence_structure.csv"), _hdr(rows), rows)
    return (; rows, cells, Bs)
end

# -----------------------------------------------------------------------------
# Study 5: commensurability control
# -----------------------------------------------------------------------------

function study_commensurability(c::Ctx, out)
    s = _sub(c.run, "commensurability")
    cells = get(s, "cells", [[12, 12, 1], [16, 16, 1], [24, 24, 1]])
    nr = Int(get(s, "n_realizations", 8))
    Bs = Float64.(get(s, "fields_T", [0.5, 1.0, 3.0, 5.0, 7.0]))
    init = String(get(s, "initialization", "field_polarized"))

    println("\n================ 5. commensurability with 3-sublattice order ================")
    println("  The 120-degree state needs L divisible by 3. Cells that are not are a")
    println("  deliberate control: a large shift signals boundary-induced frustration.")

    rows = NamedTuple[]
    for cell in cells
        label = SV.sv_cell_label(cell)
        L = SV.sv_tuple3(cell)[1]
        t = @elapsed for r in 0:(nr - 1)
            sw = _sweep(c, cell, r, init, Bs)
            for (i, B) in enumerate(Bs)
                push!(rows, (; cell=label, L, divisible_by_3=(L % 3 == 0),
                    realization=r, field_T=B, M=sw.M[i], E_per_site=sw.E[i]))
            end
        end
        @printf("  %-9s L%%3=%d  %6.1f s\n", label, L % 3, t)
    end

    println()
    @printf("  %8s %-9s %6s %12s %12s %12s\n", "B(T)", "cell", "L%3", "mean M", "sem", "E/site")
    for B in Bs
        for cell in cells
            label = SV.sv_cell_label(cell)
            sel = [x for x in rows if x.cell == label && x.field_T == B]
            isempty(sel) && continue
            d = _dist([x.M for x in sel])
            @printf("  %8.2f %-9s %6d %12.5f %12.5f %12.6f\n", B, label,
                    SV.sv_tuple3(cell)[1] % 3, d.mean, d.sem,
                    mean(x.E_per_site for x in sel))
        end
        println()
    end

    _write_csv(joinpath(out, "mvh_convergence_commensurability.csv"), _hdr(rows), rows)
    return (; rows, cells, Bs)
end

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------

function make_figure(path, sc, ini, hys, str, com)
    mkpath(dirname(path))
    fig = Figure(size=(1500, 1000))
    pal = [:royalblue, :darkorange, :seagreen, :orchid, :grey40]

    # A: scatter vs system size, log-log, against the 1/sqrt(N) reference.
    ax = Axis(fig[1, 1]; xlabel="sites", ylabel="std of M over realizations (μB/site)",
        xscale=log10, yscale=log10, title="A. Realization scatter vs cell size")
    if sc !== nothing
        for (k, B) in enumerate(sc.Bs)
            ns, sd = Float64[], Float64[]
            for cell in sc.cells
                st = only([x for x in sc.stats
                           if x.cell == SV.sv_cell_label(cell) && x.field_T == B])
                st.std > 0 || continue
                push!(ns, Float64(st.nsites)); push!(sd, st.std)
            end
            length(ns) >= 2 || continue
            scatterlines!(ax, ns, sd; color=pal[mod1(k, 5)], label=@sprintf("%.1f T", B))
        end
        n0 = Float64(prod(SV.sv_tuple3(sc.cells[1])))
        ns = range(n0, Float64(prod(SV.sv_tuple3(sc.cells[end]))); length=20)
        ref = 0.03 .* sqrt.(n0 ./ ns)
        lines!(ax, collect(ns), ref; color=:black, linestyle=:dash, label="1/√N reference")
        axislegend(ax; position=:lb, labelsize=9, framevisible=false)
    end

    # B: initialization comparison.
    ax = Axis(fig[1, 2]; xlabel="B (T)", ylabel="M (μB/site)",
        title="B. Initialization sensitivity (largest cell)")
    if ini !== nothing
        big = SV.sv_cell_label(ini.cells[end])
        for (k, init) in enumerate(ini.inits)
            ms = [mean([x.M for x in ini.rows
                        if x.cell == big && x.init == init && x.field_T == B]) for B in ini.Bs]
            scatterlines!(ax, ini.Bs, ms; color=pal[mod1(k, 5)], label=init)
        end
        axislegend(ax; position=:rb, labelsize=9, framevisible=false)
    end

    # C: hysteresis loop.
    ax = Axis(fig[2, 1]; xlabel="B (T)", ylabel="M_up − M_down (μB/site)",
        title="C. Field-history dependence")
    hlines!(ax, [0.0]; color=:black, linestyle=:dash)
    if hys !== nothing
        for (k, cell) in enumerate(hys.cells)
            label = SV.sv_cell_label(cell)
            sel = [x for x in hys.summ if x.cell == label]
            scatterlines!(ax, [x.field_T for x in sel], [x.M_up - x.M_down for x in sel];
                color=pal[mod1(k, 5)], label=label)
        end
        axislegend(ax; position=:rb, labelsize=9, framevisible=false)
    end

    # D: structure-factor width. Flat width_rlu => xi resolved.
    ax = Axis(fig[2, 2]; xlabel="1/L", ylabel="structure factor width (rlu)",
        title="D. Texture width vs cell size (flat ⇒ ξ resolved)")
    if str !== nothing
        for (k, B) in enumerate(str.Bs)
            xs, ws = Float64[], Float64[]
            for cell in str.cells
                label = SV.sv_cell_label(cell)
                sel = [x for x in str.rows if x.cell == label && x.field_T == B]
                isempty(sel) && continue
                push!(xs, mean(x.inv_L for x in sel))
                push!(ws, mean(x.width_rlu for x in sel))
            end
            length(xs) >= 2 || continue
            scatterlines!(ax, xs, ws; color=pal[mod1(k, 5)], label=@sprintf("%.1f T", B))
        end
        lines!(ax, [0.0, 0.09], [0.0, 0.09]; color=:black, linestyle=:dash,
            label="width = 1/L (resolution limit)")
        axislegend(ax; position=:lt, labelsize=9, framevisible=false)
    end

    # E: commensurability control. L not divisible by 3 cannot host the
    # 120-degree state cleanly, so it is drawn dashed.
    ax = Axis(fig[3, 1:2]; xlabel="B (T)", ylabel="M (μB/site)", xticks=0:1:7,
        title="E. Commensurability control — dashed cells are NOT divisible by 3")
    if com !== nothing
        for (k, cell) in enumerate(com.cells)
            label = SV.sv_cell_label(cell)
            L = SV.sv_tuple3(cell)[1]
            ms, es = Float64[], Float64[]
            for B in com.Bs
                sel = [x.M for x in com.rows if x.cell == label && x.field_T == B]
                isempty(sel) && continue
                d = _dist(sel)
                push!(ms, d.mean); push!(es, d.sem)
            end
            length(ms) == length(com.Bs) || continue
            errorbars!(ax, com.Bs, ms, es; color=pal[mod1(k, 5)], whiskerwidth=8)
            scatterlines!(ax, com.Bs, ms; color=pal[mod1(k, 5)],
                linestyle=(L % 3 == 0 ? :solid : :dash),
                label=@sprintf("%s (L%%3=%d)", label, L % 3))
        end
        axislegend(ax; position=:rb, labelsize=9, framevisible=false)
    end

    Label(fig[0, :], "M(H) convergence and protocol diagnostics (T = 0, minimize_energy!)";
        fontsize=16, font=:bold)
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
    params, overrides = _apply_param_overrides(params, run)

    @info "M(H) convergence diagnostics" sunny_version=SV.sv_try_pkgversion(Sunny) params_path=path diag_path base
    if !isempty(overrides)
        println("\n[run.param_overrides] applied (NOT the canonical best fit):")
        for s in overrides
            println("  ", s)
        end
    end
    c = _ctx(params, controls, run)
    @printf("\nJ1=%.6g J2=%.6g gzz=%.6g sigma_J=%.6g sigma_gzz=%.6g   T=%.3g K (kT=%.5g meV)\n",
            params.J1_meV, params.J2_meV, params.gzz, params.sigma_J, params.sigma_gzz,
            Float64(get(run, "temperature_K", 0.42)), c.kT)
    @printf("All studies at T=0 via minimize_energy!(maxiters=%d), seed_dims=%s, moment_sign=%.0f\n",
            c.maxiters, SV.sv_cell_label(c.seed_dims), c.moment_sign)

    out_tab = SV.sv_repo_path(REPO_ROOT, controls["paths"]["table_subdir"])
    out_fig = SV.sv_repo_path(REPO_ROOT, controls["paths"]["figure_subdir"])
    mkpath(out_tab); mkpath(out_fig)

    t_all = @elapsed begin
        sc  = _on(run, "scaling")          ? study_scaling(c, out_tab)          : nothing
        ini = _on(run, "initialization")   ? study_initialization(c, out_tab)   : nothing
        hys = _on(run, "hysteresis")       ? study_hysteresis(c, out_tab)       : nothing
        str = _on(run, "structure")        ? study_structure(c, out_tab)        : nothing
        com = _on(run, "commensurability") ? study_commensurability(c, out_tab) : nothing
    end

    fig_path = make_figure(joinpath(out_fig, "mvh_convergence_diagnostics.png"),
                           sc, ini, hys, str, com)

    println("\n================ recommendation ================")
    if ini !== nothing
        if ini.worst < 0.005
            @printf("  Protocol: initialization independent to %.4f uB — field_polarized +\n", ini.worst)
            println("  minimize_energy! is safe as the optimizer inner loop.")
        else
            @printf("  Protocol: initialization changes M by up to %.4f uB. Use the\n", ini.worst)
            println("  lowest-energy protocol reported above, and re-check the optimum.")
        end
    end
    if hys !== nothing
        wp = maximum(values(hys.worst_nonzero))
        if wp < 0.005
            @printf("  History: reversible for B > 0 (worst %.4f uB). The B = 0 gap is a\n", wp)
            println("  symmetry artifact and irrelevant to a fit over the measured range.")
        else
            @printf("  History: hysteretic at the %.4f uB level for B > 0 — a real bias that\n", wp)
            println("  realization averaging will NOT remove. Prefer the larger cell.")
        end
    end
    if sc !== nothing
        heavy = any((isfinite(x.skew) && abs(x.skew) > 1.0) ||
                    (isfinite(x.excess_kurtosis) && x.excess_kurtosis > 1.5) for x in sc.stats)
        println(heavy ?
            "  Statistics: heavy-tailed realization distribution — rare regions matter,\n  so use many more realizations and quote quantiles, not mean +- std." :
            "  Statistics: realization distribution is consistent with Gaussian (and in fact\n  slightly light-tailed), so mean +- sem over a modest realization count is fair.\n  Note 16 realizations only constrains skew to about +-0.6, so this is a weak test.")
        # Empirical scatter exponent: sigma ~ N^(-alpha). alpha = 0.5 is pure
        # self-averaging. alpha < 0.5 means realization averaging buys more per
        # unit compute than growing the cell.
        cs = sc.cells
        if length(cs) >= 2
            n1, n2 = prod(SV.sv_tuple3(cs[1])), prod(SV.sv_tuple3(cs[end]))
            als = Float64[]
            for B in sc.Bs
                s1 = only([x for x in sc.stats if x.cell == SV.sv_cell_label(cs[1]) && x.field_T == B]).std
                s2 = only([x for x in sc.stats if x.cell == SV.sv_cell_label(cs[end]) && x.field_T == B]).std
                (s1 > 0 && s2 > 0) && push!(als, log(s1 / s2) / log(n2 / n1))
            end
            if !isempty(als)
                a = mean(als)
                @printf("  Scaling: scatter ~ N^-%.2f (pure self-averaging would be N^-0.50).\n", a)
                if a < 0.45
                    println("  Since scatter falls SLOWER than 1/sqrt(N) while cost grows faster than")
                    println("  linearly, averaging more small cells beats growing the cell. Bound the")
                    println("  cell size by the hysteresis/texture bias, then buy accuracy with")
                    println("  realizations - which also parallelize trivially.")
                end
            end
        end
    end
    if str !== nothing
        ws = [mean(x.width_rlu for x in str.rows if x.cell == SV.sv_cell_label(cl))
              for cl in str.cells]
        Ls = [SV.sv_tuple3(cl)[1] for cl in str.cells]
        flat = length(ws) >= 2 && (maximum(ws) - minimum(ws)) < 0.2 * mean(ws)
        resolution_limited = all(i -> ws[i] < 3 / Ls[i], eachindex(ws))
        if flat && !resolution_limited
            @printf("  Texture: structure-factor width is L independent at %.3f rlu, far above\n", mean(ws))
            println("  the 1/L resolution floor. The correlation length is RESOLVED and short")
            println("  (order one lattice constant), so the cell is not limiting the texture.")
        elseif resolution_limited
            println("  Texture: width tracks the 1/L resolution floor — the correlation length")
            println("  EXCEEDS the cell. Enlarge the cell before trusting these results.")
        end
    end
    @printf("\nTotal compute: %.1f s\n", t_all)
    println("\nWrote tables to: ", out_tab)
    println("Wrote figure to: ", fig_path)
    return (; sc, ini, hys, str, com, fig_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
