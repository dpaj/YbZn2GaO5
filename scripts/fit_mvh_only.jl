#!/usr/bin/env julia

# Fit M(H) alone, and test whether the fit needs a Van Vleck term.
#
# Optimizes only what the landscape map showed M(H) can constrain
# (scripts/map_mvh_landscape.jl): J1, gzz, sigma_gzz. sigma_J is held fixed
# because M(H) is flat in it over its entire physical range, so letting an
# optimizer wander along that direction would be meaningless. J2 is held fixed
# too, being weakly constrained and preferring zero.
#
# Van Vleck test. The model is linear in both the moment amplitude A_M and the
# Van Vleck slope chi_vv, so both are profiled out by nonnegative least squares at
# every evaluation. That makes the comparison clean: refit the same Sunny curves
# with chi_vv forced to zero (a one-parameter amplitude fit) and see whether the
# rms degrades by more than the reproducibility floor. If it does not, the term is
# not needed; if it does, it is.
#
# Run WITH THREADS:
#   julia -t auto --project=. scripts/fit_mvh_only.jl

using Printf
using Statistics
using LinearAlgebra
using Optim
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

function _deepmerge!(dest::Dict, src::Dict)
    for (k, v) in src
        if v isa Dict && haskey(dest, k) && dest[k] isa Dict
            _deepmerge!(dest[k], v)
        else
            dest[k] = v
        end
    end
    return dest
end

_repo_path(root, p) = isabspath(p) ? normpath(p) : normpath(joinpath(root, splitpath(p)...))

function _load_controls(repo_root)
    rel = isempty(ARGS) || isempty(ARGS[1]) ?
        get(ENV, "SUNNY_MVH_FIT_CONTROLS", "configs/mvh_fit_controls.toml") : ARGS[1]
    p = _repo_path(repo_root, rel)
    isfile(p) || error("Could not find controls: $p")
    diag = load_toml_config(p)
    base = _repo_path(repo_root, get(get(diag, "paths", Dict{String,Any}()),
        "base_controls_toml", "configs/sunny_validation_controls.toml"))
    controls = load_toml_config(base)
    haskey(diag, "control_overrides") && _deepmerge!(controls, diag["control_overrides"])
    if haskey(diag, "paths")
        for k in ("figure_subdir", "table_subdir")
            haskey(diag["paths"], k) && (controls["paths"][k] = diag["paths"][k])
        end
    end
    return (; diag, controls, diag_path=p, base)
end

function _apply_overrides(params, run::Dict)
    ov = get(run, "param_overrides", Dict{String,Any}())
    (ov isa Dict && !isempty(ov)) || return params
    kv = Dict{Symbol,Float64}()
    for k in keys(ov)
        sym = Symbol(k)
        hasproperty(params, sym) || error("'$k' is not a model parameter")
        kv[sym] = Float64(ov[k])
    end
    return merge(params, NamedTuple(kv))
end

_setp(params, kv::Dict{Symbol,Float64}) = merge(params, NamedTuple(kv))

# -----------------------------------------------------------------------------
# Scale fits: with and without a Van Vleck term, from the SAME Sunny curve.
# -----------------------------------------------------------------------------

"One-parameter nonnegative amplitude fit, i.e. chi_vv forced to zero."
function _amplitude_only(y, x)
    m = isfinite.(y) .& isfinite.(x)
    any(m) || return 0.0
    s = sum(x[m] .^ 2)
    return s > 0 ? max(0.0, sum(y[m] .* x[m]) / s) : 0.0
end

_rms(r) = (ok = isfinite.(r); count(ok) > 0 ? sqrt(mean(r[ok] .^ 2)) : NaN)

"""
Both scale conventions for one raw Sunny curve.
`free`  : A_M and chi_vv fitted jointly (nonnegative least squares)
`novv`  : chi_vv forced to 0, amplitude only
"""
function _both_fits(raw, Bs, M_exp)
    a, b = SV.sv_best_two_component_scale(M_exp, raw, Bs)
    model_free = a .* raw .+ b .* Bs
    a0 = _amplitude_only(M_exp, raw)
    model_novv = a0 .* raw
    return (; A_M=a, chi_vv=(a > 0 ? b / a : NaN), model_free,
              rms_free=_rms(model_free .- M_exp),
              A_M_novv=a0, model_novv, rms_novv=_rms(model_novv .- M_exp),
              vv_at_max=b * maximum(Bs))
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

function main()
    (; diag, controls, diag_path, base) = _load_controls(REPO_ROOT)
    run = get(diag, "run", Dict{String,Any}())
    (; params, path) = SV.sv_load_params(REPO_ROOT, controls)
    p_byeye = _apply_overrides(params, run)

    nB = Int(get(run, "nB", 18))
    Bs = collect(range(Float64(get(run, "Bmin_T", 0.2)), Float64(get(run, "Bmax_T", 6.8)); length=nB))
    M_exp = SV.sv_mvh_target(REPO_ROOT, controls, Bs)
    nR = Int(get(run, "n_realizations", 16))
    cell = get(run, "cell_size", [12, 12, 1])
    seed_dims = get(run, "seed_dims", [3, 3, 1])
    maxiters = Int(get(run, "minimize_maxiters", 30_000))

    @info "M(H)-only fit and Van Vleck test" sunny_version=SV.sv_try_pkgversion(Sunny) params_path=path diag_path base
    @printf("\nThreads %d.  Cell %s, %d realizations, %d fields over %.3g-%.3g T.\n",
            Threads.nthreads(), SV.sv_cell_label(cell), nR, nB, first(Bs), last(Bs))
    Threads.nthreads() == 1 && @warn "Single-threaded; relaunch with `julia -t auto`."

    curve(p; realizations=0:(nR - 1), c=cell) = SV.sv_mvh_curve(p, controls;
        cell_size=c, seed_dims, realizations, Bs, maxiters, threaded=true).M_mean

    # Reproducibility floor: swap the realization set at fixed parameters.
    r_a = curve(p_byeye; realizations=0:(nR - 1))
    r_b = curve(p_byeye; realizations=nR:(2nR - 1))
    floor = abs(_both_fits(r_a, Bs, M_exp).rms_free - _both_fits(r_b, Bs, M_exp).rms_free)
    @printf("Reproducibility floor (realization set swap): %.5f uB\n", floor)

    # ---- optimize only the constrained parameters ---------------------------
    freep = Symbol.(String.(get(run, "free_parameters", ["J1_meV", "gzz", "sigma_gzz"])))
    x0 = [Float64(getproperty(p_byeye, s)) for s in freep]
    lo = Float64.(get(run, "lower_bounds", [0.10, 2.0, 0.0]))
    hi = Float64.(get(run, "upper_bounds", [0.45, 6.0, 2.0]))
    @printf("\nOptimizing %s (sigma_J = %.4g and J2 = %.4g held fixed: M(H) cannot constrain them)\n",
            join(String.(freep), ", "), p_byeye.sigma_J, p_byeye.J2_meV)

    nev = Ref(0)
    function obj(x)
        xc = clamp.(x, lo, hi)
        pen = 1e3 * sum(abs, xc .- x)          # soft box constraint
        p = _setp(p_byeye, Dict(zip(freep, xc)))
        f = _both_fits(curve(p), Bs, M_exp)
        nev[] += 1
        nev[] % 10 == 0 && @printf("    eval %3d: %s -> rms %.5f\n", nev[],
            join([@sprintf("%s=%.4g", s, v) for (s, v) in zip(freep, xc)], " "), f.rms_free)
        return (isfinite(f.rms_free) ? f.rms_free : 1e3) + pen
    end

    t_opt = @elapsed res = Optim.optimize(obj, x0, NelderMead(),
        Optim.Options(iterations=Int(get(run, "max_iterations", 120)),
                      g_tol=1e-8, show_trace=false))
    xbest = clamp.(Optim.minimizer(res), lo, hi)
    p_best = _setp(p_byeye, Dict(zip(freep, xbest)))
    @printf("\n%d evaluations in %.0f s. Optimizer: %s\n", nev[], t_opt,
            Optim.converged(res) ? "converged" : "hit iteration limit")
    for (s, v0, v) in zip(freep, x0, xbest)
        @printf("  %-11s %8.4g -> %8.4g\n", String(s), v0, v)
    end

    # ---- evaluate both parameter sets, both scale conventions ---------------
    raw_byeye = r_a
    raw_best = curve(p_best)
    f_byeye = _both_fits(raw_byeye, Bs, M_exp)
    f_best = _both_fits(raw_best, Bs, M_exp)

    muB = 0.05788381806
    bsat(p) = 0.5 * max(0.0, 9 * p.J1_meV, 8 * (p.J1_meV + p.J2_meV)) / (p.gzz * muB)

    println("\n================ Van Vleck test ================")
    println("  Same Sunny curves, refit with chi_vv free and with chi_vv forced to zero.")
    @printf("  %-22s %10s %10s %10s %10s %12s\n",
            "parameter set", "rms free", "rms novv", "d rms", "chi_vv", "VV at 6.8T")
    for (nm, f) in (("by-eye neutron", f_byeye), ("M(H)-optimized", f_best))
        @printf("  %-22s %10.5f %10.5f %10.5f %10.4f %12.4f\n",
                nm, f.rms_free, f.rms_novv, f.rms_novv - f.rms_free, f.chi_vv, f.vv_at_max)
    end
    println()
    for (nm, f) in (("by-eye neutron", f_byeye), ("M(H)-optimized", f_best))
        d = f.rms_novv - f.rms_free
        @printf("  %s: dropping Van Vleck costs %.5f uB = %.1f x the floor -> %s\n",
                nm, d, d / max(floor, 1e-12),
                d < floor ? "NOT needed" : d < 2 * floor ? "marginal" : "NEEDED")
    end

    println("\n================ parameter comparison ================")
    @printf("  %-12s %10s %10s\n", "", "by-eye", "M(H) fit")
    for s in (:J1_meV, :J2_meV, :sigma_J, :gzz, :sigma_gzz)
        @printf("  %-12s %10.4g %10.4g\n", String(s),
                Float64(getproperty(p_byeye, s)), Float64(getproperty(p_best, s)))
    end
    @printf("  %-12s %10.2f %10.2f\n", "B_sat (T)", bsat(p_byeye), bsat(p_best))
    @printf("  %-12s %10.4f %10.4f\n", "A_M", f_byeye.A_M, f_best.A_M)
    @printf("  %-12s %10.4f %10.4f\n", "chi_vv", f_byeye.chi_vv, f_best.chi_vv)
    @printf("  %-12s %10.5f %10.5f\n", "rms", f_byeye.rms_free, f_best.rms_free)

    # ---- confirm the optimum at a larger cell with different seeds ----------
    val = nothing
    if Bool(get(run, "validate_large_cell", true))
        vc = get(run, "validation_cell_size", [36, 36, 1])
        vr = Int(get(run, "validation_realization_offset", 100))
        vn = Int(get(run, "validation_n_realizations", 4))
        @printf("\nValidating at %s with realizations %d:%d (different seeds)...\n",
                SV.sv_cell_label(vc), vr, vr + vn - 1)
        rv = curve(p_best; realizations=vr:(vr + vn - 1), c=vc)
        fv = _both_fits(rv, Bs, M_exp)
        @printf("  %s: rms free %.5f (novv %.5f), A_M %.4f, chi_vv %.4f\n",
                SV.sv_cell_label(vc), fv.rms_free, fv.rms_novv, fv.A_M, fv.chi_vv)
        val = (; cell=SV.sv_cell_label(vc), raw=rv, f=fv)
    end

    # ---- outputs -----------------------------------------------------------
    out_tab = SV.sv_repo_path(REPO_ROOT, controls["paths"]["table_subdir"])
    out_fig = SV.sv_repo_path(REPO_ROOT, controls["paths"]["figure_subdir"])
    mkpath(out_tab); mkpath(out_fig)

    csv = joinpath(out_tab, "mvh_fit_comparison.csv")
    open(csv, "w") do io
        println(io, "B_T,M_exp,raw_byeye,model_byeye_free,model_byeye_novv,raw_best,model_best_free,model_best_novv,vv_byeye,vv_best")
        bb = f_byeye.chi_vv * f_byeye.A_M
        bs = f_best.chi_vv * f_best.A_M
        for i in eachindex(Bs)
            println(io, join([@sprintf("%.10g", v) for v in (
                Bs[i], M_exp[i], raw_byeye[i], f_byeye.model_free[i], f_byeye.model_novv[i],
                raw_best[i], f_best.model_free[i], f_best.model_novv[i],
                bb * Bs[i], bs * Bs[i])], ","))
        end
    end

    summ = joinpath(out_tab, "mvh_fit_summary.csv")
    open(summ, "w") do io
        println(io, "set,J1_meV,J2_meV,sigma_J,gzz,sigma_gzz,B_sat_T,A_M,chi_vv,rms_free,rms_novv,d_rms,floor,cell,n_realizations")
        for (nm, p, f, c, n) in (("by-eye", p_byeye, f_byeye, SV.sv_cell_label(cell), nR),
                                 ("mvh_fit", p_best, f_best, SV.sv_cell_label(cell), nR))
            println(io, join(Any[nm, p.J1_meV, p.J2_meV, p.sigma_J, p.gzz, p.sigma_gzz,
                bsat(p), f.A_M, f.chi_vv, f.rms_free, f.rms_novv,
                f.rms_novv - f.rms_free, floor, c, n], ","))
        end
        if val !== nothing
            f = val.f
            println(io, join(Any["mvh_fit_validated", p_best.J1_meV, p_best.J2_meV, p_best.sigma_J,
                p_best.gzz, p_best.sigma_gzz, bsat(p_best), f.A_M, f.chi_vv, f.rms_free,
                f.rms_novv, f.rms_novv - f.rms_free, floor, val.cell,
                Int(get(run, "validation_n_realizations", 4))], ","))
        end
    end

    # ---- figure ------------------------------------------------------------
    data = SV.sv_read_magnetization_csv(SV.sv_repo_path(REPO_ROOT, controls["paths"]["magnetization_csv"]))
    fig = Figure(size=(1450, 900))

    ax = Axis(fig[1, 1]; xlabel="B (T)", ylabel="M (μB / Yb)", xticks=0:1:7,
        title="A. M(H): experiment vs by-eye neutron parameters vs M(H)-only fit")
    scatter!(ax, data.B_T, data.M_muB_per_Yb; markersize=3, color=:black, label="experiment (digitized)")
    lines!(ax, Bs, f_byeye.model_free; color=:darkorange, linewidth=2.5,
        label=@sprintf("by-eye neutron (rms %.4f)", f_byeye.rms_free))
    lines!(ax, Bs, f_best.model_free; color=:royalblue, linewidth=2.5,
        label=@sprintf("M(H)-only fit (rms %.4f)", f_best.rms_free))
    axislegend(ax; position=:rb, framevisible=false, labelsize=10)

    ax = Axis(fig[1, 2]; xlabel="B (T)", ylabel="model − experiment (μB / Yb)", xticks=0:1:7,
        title="B. Residuals")
    hlines!(ax, [0.0]; color=:black, linestyle=:dash)
    band!(ax, Bs, fill(-floor, length(Bs)), fill(floor, length(Bs));
        color=(:grey, 0.25), label="reproducibility floor")
    lines!(ax, Bs, f_byeye.model_free .- M_exp; color=:darkorange, linewidth=2.5, label="by-eye")
    lines!(ax, Bs, f_best.model_free .- M_exp; color=:royalblue, linewidth=2.5, label="M(H) fit")
    axislegend(ax; position=:lb, framevisible=false, labelsize=10)

    ax = Axis(fig[2, 1]; xlabel="B (T)", ylabel="M (μB / Yb)", xticks=0:1:7,
        title="C. Van Vleck contribution and the model without it")
    scatter!(ax, data.B_T, data.M_muB_per_Yb; markersize=3, color=:black, label="experiment")
    lines!(ax, Bs, f_best.model_free; color=:royalblue, linewidth=2.5, label="M(H) fit, χ_vv free")
    lines!(ax, Bs, f_best.model_novv; color=:royalblue, linestyle=:dash, linewidth=2,
        label=@sprintf("M(H) fit, χ_vv ≡ 0 (rms %.4f)", f_best.rms_novv))
    lines!(ax, Bs, f_best.A_M * f_best.chi_vv .* Bs; color=:crimson, linestyle=:dot, linewidth=2,
        label=@sprintf("its χ_vv part (%.3f μB at %.1f T)", f_best.A_M * f_best.chi_vv * last(Bs), last(Bs)))
    lines!(ax, Bs, f_byeye.A_M * f_byeye.chi_vv .* Bs; color=:darkorange, linestyle=:dot, linewidth=2,
        label=@sprintf("by-eye χ_vv part (%.3f μB)", f_byeye.A_M * f_byeye.chi_vv * last(Bs)))
    axislegend(ax; position=:lt, framevisible=false, labelsize=9)

    ax = Axis(fig[2, 2]; xlabel="B (T)", ylabel="model − experiment (μB / Yb)", xticks=0:1:7,
        title="D. Cost of forcing χ_vv = 0")
    hlines!(ax, [0.0]; color=:black, linestyle=:dash)
    band!(ax, Bs, fill(-floor, length(Bs)), fill(floor, length(Bs)); color=(:grey, 0.25))
    lines!(ax, Bs, f_best.model_free .- M_exp; color=:royalblue, linewidth=2.5, label="χ_vv free")
    lines!(ax, Bs, f_best.model_novv .- M_exp; color=:royalblue, linestyle=:dash, linewidth=2,
        label="χ_vv ≡ 0")
    lines!(ax, Bs, f_byeye.model_novv .- M_exp; color=:darkorange, linestyle=:dash, linewidth=2,
        label="by-eye, χ_vv ≡ 0")
    axislegend(ax; position=:lb, framevisible=false, labelsize=9)

    Label(fig[0, :], @sprintf("M(H)-only fit vs by-eye neutron parameters — %s cell, %d realizations, A_M profiled out (shape only)",
            SV.sv_cell_label(cell), nR); fontsize=15, font=:bold)
    figp = joinpath(out_fig, "mvh_fit_comparison.png")
    save(figp, fig)

    println("\nWrote:\n  ", csv, "\n  ", summ, "\n  ", figp)
    return (; p_byeye, p_best, f_byeye, f_best, floor, val, figp)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
