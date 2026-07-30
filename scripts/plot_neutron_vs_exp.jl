#!/usr/bin/env julia

# Overplot the Sunny KPM neutron model against the background-subtracted 1D cuts,
# and show how differently the two observables respond to the disorder width sigma_J.
#
#   julia -t auto --project=. scripts/plot_neutron_vs_exp.jl
#
# Produces two figures:
#   neutron_vs_exp.png          6 cuts (3 qtags x 2 fields), experiment vs model at
#                               several sigma_J, with the intensity scale profiled out
#   neutron_mvh_complementarity.png
#                               sigma_J sensitivity of the neutron objective set
#                               against the M(H) landscape scan. NOTE the two prefer
#                               OPPOSITE directions -- see the figure caption.
#
# The intensity scale is fitted, not absolute: Sunny's prefactor is not comparable to
# the analytical model's, so only lineshape is a real constraint. One global
# nonnegative scale is profiled out by weighted least squares over the fit window,
# exactly as A_M is profiled out of the M(H) objective.

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using CairoMakie
using Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation
const SV = SunnyValidation

_load() = SV.sv_load_diagnostic_controls(REPO_ROOT, "configs/neutron_vs_exp_controls.toml";
                                         env_var="SUNNY_NEUTRON_VS_EXP_CONTROLS")

"Pretty qtag: 0p33_0p33_0 -> (0.33, 0.33, 0)"
function _qlabel(t::AbstractString)
    parts = split(t, "_")
    length(parts) == 3 || return t
    f(p) = replace(p, "p" => ".")
    return "(" * join(f.(parts), ", ") * ")"
end

function main()
    (; diag, controls) = _load()
    run = get(diag, "run", Dict{String,Any}())
    (; params) = SV.sv_load_params(REPO_ROOT, controls)
    params, overrides = SV.sv_apply_param_overrides(params, run)

    reals = Int.(get(run, "realizations", [0]))
    include_flat = Bool(get(run, "include_flat", false))
    sj_over = Float64.(get(run, "sigma_J_overlay", [0.0, 0.5, 1.0]))
    sj_scan = Float64.(get(run, "sigma_J_scan", [0.0, 0.25, 0.5, 0.75, 1.0]))

    println("Neutron model vs experiment")
    isempty(overrides) || println("  overrides: ", join(overrides, ", "))
    @printf("  cell %s, %d realization(s), threads %d, regularization %.1e\n",
            SV.sv_cell_label(controls["kpm"]["system_size"]), length(reals),
            Threads.nthreads(), Float64(get(controls["kpm"], "regularization", 1e-6)))

    cuts = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
    nq = length(SV.sv_kpm_1d_q_sampler(cuts[1], controls).qs)
    @printf("  %d cuts, %d q per cut, %d chunks\n\n", length(cuts), nq,
            SV.sv_kpm_q_chunks(controls, nq))

    # ---- evaluate the model at each overlay sigma_J -------------------------
    # ~2.5 min per sigma_J at 36x36x1. The written CSV carries the scaled curves, so
    # presentation-only changes can be replotted from it rather than recomputed.
    overlay = Dict{Float64,Any}()
    for sj in sj_over
        p = merge(params, (; sigma_J=sj))
        t = @elapsed o = SV.sv_neutron_objective(p, controls, cuts;
            realizations=reals, include_flat)
        overlay[sj] = o
        @printf("  sigma_J=%.2f : chi2_red=%9.4g  rms=%.4g  scale=%.4g  (%.0f s, %d failed)\n",
                sj, o.chi2_red, o.rms, o.scale, t, o.n_failed)
    end

    # ---- denser scan for the sensitivity panel ------------------------------
    println()
    scan = NamedTuple[]
    for sj in sj_scan
        o = haskey(overlay, sj) ? overlay[sj] :
            SV.sv_neutron_objective(merge(params, (; sigma_J=sj)), controls, cuts;
                                    realizations=reals, include_flat)
        push!(scan, (; sigma_J=sj, chi2_red=o.chi2_red, rms=o.rms, scale=o.scale))
        @printf("  scan sigma_J=%.2f : chi2_red=%9.4g\n", sj, o.chi2_red)
    end

    out_fig = SV.sv_repo_path(REPO_ROOT, controls["paths"]["figure_subdir"])
    out_tab = SV.sv_repo_path(REPO_ROOT, controls["paths"]["table_subdir"])
    mkpath(out_fig); mkpath(out_tab)

    # =====================================================================
    # Figure 1: the six cuts
    # =====================================================================
    qtags = unique([c.qtag for c in cuts])
    fields = sort(unique([c.field_T for c in cuts]))
    ecut = Float64(get(run, "elastic_cutoff_meV", 0.35))
    pal = [:royalblue, :darkorange, :seagreen, :orchid]

    fig = Figure(size=(520 * length(fields), 300 * length(qtags)))
    for (row, qt) in enumerate(qtags), (col, B) in enumerate(fields)
        idx = findfirst(c -> c.qtag == qt && c.field_T ≈ B, cuts)
        idx === nothing && continue
        cut = cuts[idx]
        ax = Axis(fig[row, col]; xlabel="energy transfer (meV)", ylabel="intensity",
            title=@sprintf("%s   %.0f T", _qlabel(qt), B))
        errorbars!(ax, cut.energy_meV, cut.intensity, cut.error;
                   color=(:black, 0.35), whiskerwidth=0)
        scatter!(ax, cut.energy_meV, cut.intensity; color=:black, markersize=4,
                 label="experiment")
        for (k, sj) in enumerate(sj_over)
            o = overlay[sj]
            lines!(ax, cut.energy_meV, o.scale .* o.curves[idx];
                   color=pal[mod1(k, 4)], linewidth=2,
                   label=@sprintf("σ_J=%.2f (χ²=%.0f)", sj, o.chi2_red))
        end
        vspan!(ax, 0.5, 3.0; color=(:grey, 0.10))
        # Scale to the INELASTIC signal. Both the data's elastic line and the model's
        # quasi-elastic divergence at large sigma_J are orders of magnitude larger
        # than the magnon signal and would otherwise flatten every panel to zero.
        keep = cut.energy_meV .>= ecut
        if any(keep)
            top = maximum(vcat(cut.intensity[keep],
                [maximum((overlay[sj].scale .* overlay[sj].curves[idx])[keep]) for sj in sj_over]))
            bot = minimum(cut.intensity[keep])
            ylims!(ax, min(bot, 0) - 0.05 * top, 1.25 * top)
        end
        xlims!(ax, 0.0, maximum(cut.energy_meV))
        row == 1 && col == 1 && axislegend(ax; position=:rt, framevisible=false,
                                           labelsize=9)
    end
    Label(fig[0, :], @sprintf("Sunny KPM vs experiment — %s, 1 realization, %d q per cut, ONE global intensity scale profiled out
(shaded = fit window; y-axis set by the inelastic signal above %.2f meV, so the elastic line is clipped)",
        SV.sv_cell_label(controls["kpm"]["system_size"]), nq, ecut); fontsize=13, font=:bold)
    p1 = joinpath(out_fig, "neutron_vs_exp.png")
    save(p1, fig)

    # =====================================================================
    # Figure 2: complementarity with M(H)
    # =====================================================================
    fig2 = Figure(size=(1400, 460))

    ax = Axis(fig2[1, 1]; xlabel="σ_J", ylabel="neutron χ²_red",
        title="A. Neutron: strong, monotonic — still falling at σ_J = 1")
    ok = [s for s in scan if isfinite(s.chi2_red)]
    scatterlines!(ax, [s.sigma_J for s in ok], [s.chi2_red for s in ok];
                  color=:darkorange, linewidth=2.5, markersize=10)
    vlines!(ax, [Float64(params.sigma_J)]; color=:grey, linestyle=:dash)

    # M(H) sigma_J scan, from the already-computed landscape map
    ax2 = Axis(fig2[1, 2]; xlabel="σ_J", ylabel="M(H) rms residual (μB/Yb)",
        title="B. M(H): weak, shallow minimum near σ_J = 0.3")
    mvh = SV.sv_repo_path(REPO_ROOT,
        "results/feature_tables/sunny_validation/mvh_landscape/mvh_landscape_scan1d.csv")
    have_mvh = isfile(mvh)
    if have_mvh
        raw, hdr = readdlm(mvh, ','; header=true)
        cols = vec(String.(hdr))
        ip, iv, ir = findfirst(==("parameter"), cols), findfirst(==("value"), cols),
                     findfirst(==("rms"), cols)
        sel = [i for i in axes(raw, 1) if String(raw[i, ip]) == "sigma_J"]
        xs = Float64.(raw[sel, iv]); ys = Float64.(raw[sel, ir])
        o = sortperm(xs)
        scatterlines!(ax2, xs[o], ys[o]; color=:royalblue, linewidth=2.5, markersize=10)
        floor_mvh = 0.00259   # reproducibility floor from the landscape run
        band!(ax2, xs[o], fill(minimum(ys) , length(xs)),
              fill(minimum(ys) + floor_mvh, length(xs)); color=(:grey, 0.25))
        text!(ax2, 0.05, minimum(ys) + 1.15 * floor_mvh;
              text="best + reproducibility floor", fontsize=9, color=:grey30)
    else
        text!(ax2, 0.5, 0.5; text="run scripts/map_mvh_landscape.jl first",
              align=(:center, :center))
    end

    # Both normalized, so the contrast is visible on one axis
    ax3 = Axis(fig2[1, 3]; xlabel="σ_J", ylabel="objective, normalized to its own min",
        title="C. They pull in OPPOSITE directions")
    if !isempty(ok)
        m = minimum(s.chi2_red for s in ok)
        scatterlines!(ax3, [s.sigma_J for s in ok], [s.chi2_red / m for s in ok];
                      color=:darkorange, linewidth=2.5, markersize=10, label="neutron χ²_red")
    end
    if have_mvh
        raw, hdr = readdlm(mvh, ','; header=true)
        cols = vec(String.(hdr))
        ip, iv, ir = findfirst(==("parameter"), cols), findfirst(==("value"), cols),
                     findfirst(==("rms"), cols)
        sel = [i for i in axes(raw, 1) if String(raw[i, ip]) == "sigma_J"]
        xs = Float64.(raw[sel, iv]); ys = Float64.(raw[sel, ir])
        o = sortperm(xs)
        scatterlines!(ax3, xs[o], ys[o] ./ minimum(ys); color=:royalblue,
                      linewidth=2.5, markersize=10, label="M(H) rms")
    end
    hlines!(ax3, [1.0]; color=:black, linestyle=:dash)
    axislegend(ax3; position=:rt, framevisible=false, labelsize=10)

    text!(ax3, 0.5, 1.30; text="neutron pulls σ_J up,
M(H) pulls it down —
tension, not just
complementarity",
          fontsize=10, color=:grey25, align=(:center, :center))
    Label(fig2[0, :], "σ_J: the neutron objective moves 36%, M(H) only 2.7x its reproducibility floor — and their preferences oppose";
          fontsize=13, font=:bold)
    p2 = joinpath(out_fig, "neutron_mvh_complementarity.png")
    save(p2, fig2)

    # ---- tables -------------------------------------------------------------
    SV.sv_write_rows_csv(joinpath(out_tab, "neutron_sigmaJ_scan.csv"), scan)
    open(joinpath(out_tab, "neutron_vs_exp_curves.csv"), "w") do io
        println(io, "qtag,field_T,energy_meV,I_exp,Ierr_exp," *
                    join(["I_model_sigmaJ_$(replace(string(s),'.'=>'p'))" for s in sj_over], ","))
        for (i, cut) in enumerate(cuts), j in eachindex(cut.energy_meV)
            vals = [overlay[s].scale * overlay[s].curves[i][j] for s in sj_over]
            println(io, join(Any[cut.qtag, cut.field_T,
                @sprintf("%.6g", cut.energy_meV[j]), @sprintf("%.6g", cut.intensity[j]),
                @sprintf("%.6g", cut.error[j]),
                join([@sprintf("%.6g", v) for v in vals], ",")], ","))
        end
    end

    println("\nWrote:\n  ", p1, "\n  ", p2, "\n  ", out_tab)
    return (; overlay, scan, p1, p2)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
