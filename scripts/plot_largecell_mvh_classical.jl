#!/usr/bin/env julia

# Plot-only companion to scripts/sunny_largecell_mvh_classical.jl.
#
# Reads the CSV that script already wrote and makes a four-panel diagnostic
# figure. No Sunny, no recomputation, so it is cheap to iterate on the plot.
#
# Run with:
#   julia --project=. scripts/plot_largecell_mvh_classical.jl
#
# Optional positional argument: the diagnostic controls TOML, so an alternative
# run directory can be plotted.

using Printf
using Statistics
using DelimitedFiles
using CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

const SV = SunnyValidation

function _repo_path(root, p)
    isabspath(p) && return normpath(p)
    return normpath(joinpath(root, splitpath(p)...))
end

function _read_csv(path)
    raw, header = readdlm(path, ','; header=true)
    cols = vec(String.(header))
    idx = Dict(c => i for (i, c) in enumerate(cols))
    return (; raw, idx)
end

_col(t, name) = t.raw[:, t.idx[name]]
_fcol(t, name) = [x isa AbstractString ? parse(Float64, x) : Float64(x) for x in _col(t, name)]
_scol(t, name) = String.(string.(_col(t, name)))

# Average a column over disorder realizations, keyed by field, for one
# (cell, sampler, variant) selection. Returns (Bs, mean, std_over_realizations).
function _by_field(t, cell, sampler, valcol; variant=nothing)
    cells, samplers = _scol(t, "cell"), _scol(t, "sampler")
    variants = haskey(t.idx, "variant") ? _scol(t, "variant") : fill("disordered", length(cells))
    Bs, vals = _fcol(t, "B_T"), _fcol(t, valcol)
    acc = Dict{Float64,Vector{Float64}}()
    for i in eachindex(Bs)
        (cells[i] == cell && samplers[i] == sampler) || continue
        (variant === nothing || variants[i] == variant) || continue
        push!(get!(acc, Bs[i], Float64[]), vals[i])
    end
    ks = sort(collect(keys(acc)))
    mu = [mean(acc[k]) for k in ks]
    sd = [length(acc[k]) > 1 ? std(acc[k]) : 0.0 for k in ks]
    return (ks, mu, sd)
end

# Scalar knob replicated on every row of a (cell, sampler, variant) group.
function _scalar(t, cell, sampler, variant, col)
    cells, samplers = _scol(t, "cell"), _scol(t, "sampler")
    variants = haskey(t.idx, "variant") ? _scol(t, "variant") : fill("disordered", length(cells))
    vals = _fcol(t, col)
    for i in eachindex(vals)
        if cells[i] == cell && samplers[i] == sampler && variants[i] == variant
            return vals[i]
        end
    end
    return NaN
end

function main()
    diag_rel = isempty(ARGS) ? "configs/sunny_largecell_mvh_classical_controls.toml" : ARGS[1]
    diag = load_toml_config(_repo_path(REPO_ROOT, diag_rel))
    controls = load_toml_config(_repo_path(REPO_ROOT,
        get(diag["paths"], "base_controls_toml", "configs/sunny_validation_controls.toml")))
    run = get(diag, "run", Dict{String,Any}())

    tab_dir = _repo_path(REPO_ROOT, diag["paths"]["table_subdir"])
    fig_dir = _repo_path(REPO_ROOT, diag["paths"]["figure_subdir"])
    csv_path = joinpath(tab_dir, "sunny_largecell_mvh_classical.csv")
    isfile(csv_path) || error("Missing $csv_path — run scripts/sunny_largecell_mvh_classical.jl first")
    t = _read_csv(csv_path)

    data = SV.sv_read_magnetization_csv(
        _repo_path(REPO_ROOT, controls["paths"]["magnetization_csv"]))

    cells = unique(_scol(t, "cell"))
    # Order cell labels by site count rather than lexically.
    nsites = Dict(_scol(t, "cell")[i] => _fcol(t, "nsites")[i] for i in eachindex(_scol(t, "cell")))
    sort!(cells; by=c -> nsites[c])
    prod_cell = cells[end]
    T_K = Float64(get(run, "temperature_K", 0.42))
    variants = haskey(t.idx, "variant") ? unique(_scol(t, "variant")) : ["disordered"]
    sort!(variants; by=v -> v == "disordered" ? 0 : 1)

    samplers = ["minimize_energy", "langevin", "langevin_then_midpoint"]
    samplers = [s for s in samplers if any(_scol(t, "sampler") .== s)]
    samp = "langevin" in samplers ? "langevin" : samplers[1]
    vcol = Dict("disordered" => :darkorange, "clean" => :royalblue)
    ovr = get(run, "param_overrides", Dict{String,Any}())
    ovr_str = isempty(ovr) ? "canonical best-fit parameters" :
        join([@sprintf("%s=%.4g", k, Float64(ovr[k])) for k in sort(collect(keys(ovr)))], ", ")

    fig = Figure(size=(1500, 980))

    # ---- Panel A: disordered vs clean, joint (A_M, chi_vv) fit --------------
    axA = Axis(fig[1, 1]; xlabel="B (T)", ylabel="M (μB / Yb)", xticks=0:1:7,
        title=@sprintf("A. %s, %s — A_M and χ_vv fitted jointly", prod_cell, samp))
    scatter!(axA, data.B_T, data.M_muB_per_Yb; markersize=3, color=:black, label="experiment")
    for v in variants
        B, mu, sd = _by_field(t, prod_cell, samp, "M_total_joint_fit_uB"; variant=v)
        isempty(B) && continue
        aj = _scalar(t, prod_cell, samp, v, "magnetization_global_scale_joint")
        cj = _scalar(t, prod_cell, samp, v, "chi_vv_joint_fit_muB_per_T")
        lines!(axA, B, mu; color=get(vcol, v, :seagreen), linewidth=2,
            label=@sprintf("%s (A_M=%.3f, χ_vv=%.3f)", v, aj, cj))
        any(sd .> 0) && band!(axA, B, mu .- sd, mu .+ sd; color=(get(vcol, v, :seagreen), 0.2))
    end
    axislegend(axA; position=:lt, framevisible=false, labelsize=10)

    # ---- Panel B: same, but chi_vv held at the analytical value -------------
    axB = Axis(fig[1, 2]; xlabel="B (T)", ylabel="M (μB / Yb)", xticks=0:1:7,
        title="B. A_M free, χ_vv held at the analytical value")
    scatter!(axB, data.B_T, data.M_muB_per_Yb; markersize=3, color=:black, label="experiment")
    for v in variants, (k, s) in enumerate(samplers)
        B, mu, _ = _by_field(t, prod_cell, s, "M_total_free_scale_uB"; variant=v)
        isempty(B) && continue
        af = _scalar(t, prod_cell, s, v, "magnetization_global_scale_free")
        lines!(axB, B, mu; color=get(vcol, v, :seagreen),
            linestyle=(:solid, :dash, :dot)[min(k, 3)], linewidth=2,
            label=@sprintf("%s / %s (A_M=%.3f)", v, s, af))
    end
    axislegend(axB; position=:lt, framevisible=false, labelsize=9)

    # ---- Panel C: residuals — the informative panel -------------------------
    axC = Axis(fig[2, 1]; xlabel="B (T)", ylabel="model − experiment (μB / Yb)", xticks=0:1:7,
        title="C. Residual: joint fit (solid) vs fixed analytical χ_vv (dashed)")
    hlines!(axC, [0.0]; color=:black, linestyle=:dash, linewidth=1)
    for v in variants
        B, mj, _ = _by_field(t, prod_cell, samp, "residual_joint_fit_uB"; variant=v)
        isempty(B) && continue
        keep = isfinite.(mj)
        lines!(axC, B[keep], mj[keep]; color=get(vcol, v, :seagreen), linewidth=2.5,
            label=@sprintf("%s joint (max |r|=%.3f)", v, maximum(abs, mj[keep])))
        # Amplitude-only fit for the same selection, to show what fixing chi_vv costs.
        _, mu, _ = _by_field(t, prod_cell, samp, "M_total_free_scale_uB"; variant=v)
        _, ex, _ = _by_field(t, prod_cell, samp, "M_exp_interp_uB_per_Yb"; variant=v)
        k2 = isfinite.(ex)
        r2 = mu[k2] .- ex[k2]
        lines!(axC, B[k2], r2; color=get(vcol, v, :seagreen), linestyle=:dash, linewidth=1.5,
            label=@sprintf("%s fixed χ_vv (max |r|=%.3f)", v, maximum(abs, r2)))
    end
    axislegend(axC; position=:lb, framevisible=false, labelsize=9)

    # ---- Panel D: components and cell-size convergence ---------------------
    axD = Axis(fig[2, 2]; xlabel="B (T)", ylabel="M (μB / Yb)", xticks=0:1:7,
        title=@sprintf("D. Components (%s, joint fit) and cell-size convergence", samp))
    scatter!(axD, data.B_T, data.M_muB_per_Yb; markersize=3, color=:black, label="experiment")
    for (k, c) in enumerate(cells)
        B, mu, _ = _by_field(t, c, samp, "M_total_joint_fit_uB"; variant="disordered")
        isempty(B) && continue
        aj = _scalar(t, c, samp, "disordered", "magnetization_global_scale_joint")
        lines!(axD, B, mu; color=[:royalblue, :darkorange, :seagreen][min(k, 3)], linewidth=2,
            label=@sprintf("%s disordered (A_M=%.3f)", c, aj))
    end
    B, vvj, _ = _by_field(t, prod_cell, samp, "M_vv_joint_fit_uB"; variant="disordered")
    aj = _scalar(t, prod_cell, samp, "disordered", "magnetization_global_scale_joint")
    _, disp, _ = _by_field(t, prod_cell, samp, "M_disp_raw_uB_per_site"; variant="disordered")
    if !isempty(B)
        lines!(axD, B, vvj; color=:crimson, linestyle=:dot, linewidth=2,
            label=@sprintf("Van Vleck part of joint fit (%.2f μB at %.0f T)", vvj[end], B[end]))
        lines!(axD, B, aj .* disp; color=:grey40, linestyle=:dash, linewidth=2,
            label="disordered-phase moment × A_M")
    end
    axislegend(axD; position=:lt, framevisible=false, labelsize=9)

    Label(fig[0, :], @sprintf("Sunny large-cell classical M(H,T) at %.2f K — minimal single disordered phase — %s", T_K, ovr_str);
        fontsize=15, font=:bold)

    mkpath(fig_dir)
    out = joinpath(fig_dir, "sunny_largecell_mvh_classical_diagnostic.png")
    save(out, fig)
    println("Wrote: ", out)
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
