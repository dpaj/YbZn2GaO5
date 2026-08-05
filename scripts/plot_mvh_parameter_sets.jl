#!/usr/bin/env julia
# M(H) against experiment for several NAMED parameter sets, including a CLEAN (no-disorder)
# variant of the same Hamiltonian.
#
#   julia -t auto --project=. scripts/plot_mvh_parameter_sets.jl
#
# The point of interest: the neutron-fitted parameters were obtained from the SPECTRA ALONE
# and never shown M(H). M(H) constrains B_sat ~ J1/gzz, and the fit lowered that ratio by
# about a third from the by-eye pair -- toward what M(H) independently prefers. If M(H)
# improves at those parameters, that is cross-observable evidence requiring NO choice of
# relative weight between the two observables, which is the part of a co-fit that cannot be
# justified by assertion.
#
# The clean curve is a control, not a claim about the published model: this Hamiltonian is
# isotropic Heisenberg with all anisotropy in the g tensor, whereas the published fits use
# XXZ with Delta ~ 1.35. It isolates what DISORDER does, nothing more.
#
# A_M and the linear (Van Vleck-like) slope are profiled out per set by nonnegative least
# squares, so panels 1-4 compare SHAPE.
#
# PANEL 5 IS NEW AND IS THE POINT OF THE UPDATE. The old comment here said the absolute scale was
# not trustworthy because the experimental normalisation was suspect. That is no longer true: the
# MPMS3 centring correction fixed it, and three independent numbers now agree at ~7 T to four
# digits -- MPMS3 0.42 K free-centring 1.6512, DynaCool 2.5 K on a DIFFERENT crystal 1.6522, and the
# Bag et al. Supplement ~1.65. So the absolute comparison is available for the first time, and it
# tests something the shape fit CANNOT: A_M being profiled out means a level error is invisible by
# construction. Panel 5 shows the model with NO A_M against the data, and the ratio versus field.

using Printf, Statistics, LinearAlgebra, CairoMakie, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const LOADED = SV.sv_load_diagnostic_controls(REPO_ROOT,
    "configs/mvh_parameter_sets_controls.toml"; env_var="YZGO_MVH_SETS_CONTROLS")
const CFG = LOADED.diag
const RUN = get(CFG, "run", Dict{String,Any}())
const controls = LOADED.controls

params0, applied = SV.sv_apply_param_overrides(
    SV.sv_load_params(REPO_ROOT, controls).params, RUN)

const CELL = Tuple(Int.(get(RUN, "cell_size", [12, 12, 1])))
const SEEDD = Tuple(Int.(get(RUN, "seed_dims", [3, 3, 1])))
const NREAL = Int(get(RUN, "n_realizations", 16))
const MAXIT = Int(get(RUN, "minimize_maxiters", 2000))
const FDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "figure_subdir",
    "results/figures/sunny_validation/mvh_parameter_sets"))
const TDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "table_subdir",
    "results/feature_tables/sunny_validation/mvh_parameter_sets"))
mkpath(FDIR); mkpath(TDIR)

const Bs = collect(range(Float64(get(RUN, "field_min_T", 0.2)),
                         Float64(get(RUN, "field_max_T", 6.8));
                         length = Int(get(RUN, "n_fields", 34))))
const M_exp = SV.sv_mvh_target(REPO_ROOT, controls, Bs)

@printf("threads %d, cell %s, %d realizations, %d fields %.2f-%.2f T\n",
        Threads.nthreads(), string(CELL), NREAL, length(Bs), first(Bs), last(Bs))
isempty(applied) || @printf("base overrides: %s\n", join(applied, ", "))

sets = get(CFG, "sets", Any[])
isempty(sets) && error("No [[sets]] in the config.")

results = NamedTuple[]
for (i, sd) in enumerate(sets)
    name = String(get(sd, "name", "set$i"))
    over = NamedTuple(Symbol(k) => Float64(v) for (k, v) in sd if k != "name")
    p = merge(params0, over)
    t = time()
    o = SV.sv_mvh_objective(p, controls, Bs, M_exp; cell_size=CELL, seed_dims=SEEDD,
            realizations=0:(NREAL - 1), maxiters=MAXIT, threaded=true)
    @printf("\n%-18s J1=%.3f sigma_J=%.2f gzz=%.2f sigma_gzz=%.2f\n",
            name, p.J1_meV, p.sigma_J, p.gzz, p.sigma_gzz)
    @printf("    rms = %.5f uB   max|res| = %.5f   A_M = %.4g   chi_vv = %.4g   J1/gzz = %.4f  (%.0f s)\n",
            o.rms, o.max_abs, o.A_M, o.chi_vv, p.J1_meV / p.gzz, time() - t)
    # CONTROL: refit with the linear term FORCED TO ZERO. The two-component fit is free to
    # buy rms with a large chi_vv, and the crystal field only allows 0.0171 +- 0.0007 uB/T,
    # so an improvement that survives only WITH a big linear term is not an improvement in
    # the spin physics. One-component nonnegative least squares is a = (raw.M)/(raw.raw).
    a1 = max(0.0, dot(o.raw, M_exp) / max(eps(), dot(o.raw, o.raw)))
    res1 = a1 .* o.raw .- M_exp
    ok1 = isfinite.(res1)
    rms1 = sqrt(mean(res1[ok1] .^ 2))
    @printf("    control, chi_vv forced to 0:  rms = %.5f uB  (A_M = %.4g)  %s\n",
            rms1, a1, rms1 < 1.3 * o.rms ? "linear term is NOT load-bearing" :
                        "*** the fit LEANS on the linear term ***")
    push!(results, (; name, params=p, obj=o, rms_nolin=rms1, A_M_nolin=a1,
                      model_nolin=a1 .* o.raw))
    flush(stdout)
end

open(joinpath(TDIR, "mvh_parameter_sets.csv"), "w") do io
    println(io, "set,J1_meV,sigma_J,gzz,sigma_gzz,J1_over_gzz,rms_uB,rms_uB_nolinear,max_abs_uB,A_M,chi_vv,B_T,M_exp,M_model,M_model_nolinear,residual")
    for r in results, j in eachindex(Bs)
        @printf(io, "%s,%.4f,%.4f,%.4f,%.4f,%.5f,%.6g,%.6g,%.6g,%.6g,%.6g,%.4f,%.6g,%.6g,%.6g,%.6g\n",
                r.name, r.params.J1_meV, r.params.sigma_J, r.params.gzz,
                r.params.sigma_gzz, r.params.J1_meV / r.params.gzz, r.obj.rms,
                r.rms_nolin, r.obj.max_abs, r.obj.A_M, r.obj.chi_vv, Bs[j],
                M_exp[j], r.obj.model[j], r.model_nolin[j], r.obj.residual[j])
    end
end

# Central-difference susceptibility. A magnetization PLATEAU would appear here as a clear
# dip toward zero, which is the published clean model's signature and is what the data
# should be checked against.
function deriv(x, y)
    d = similar(y)
    for i in eachindex(y)
        lo = max(firstindex(y), i - 1); hi = min(lastindex(y), i + 1)
        d[i] = (y[hi] - y[lo]) / (x[hi] - x[lo])
    end
    return d
end

cols = [:crimson, :dodgerblue, :seagreen, :darkorange, :purple]
fig = Figure(size = (1500, 1480))
Label(fig[0, 1:2],
      "YbZn2GaO5 -- M(H) at 0.42 K.  Parameters fitted to the NEUTRON SPECTRA ALONE, " *
      "tested against magnetization";
      fontsize = 17, font = :bold)

# ---- 1  M(H) ---------------------------------------------------------------------
ax1 = Axis(fig[1, 1], xlabel = "field (T)", ylabel = "M (uB / Yb)",
           title = "A_M and the linear slope profiled out per set (shape comparison)")
scatter!(ax1, Bs, M_exp; color = :black, markersize = 9, label = "experiment")
for (k, r) in enumerate(results)
    lines!(ax1, Bs, r.obj.model; color = cols[mod1(k, length(cols))], linewidth = 2.6,
           label = "$(r.name)  (rms $(round(r.obj.rms; sigdigits=3)))")
end
axislegend(ax1; position = :rb, labelsize = 10)

# ---- 2  residuals ----------------------------------------------------------------
ax2 = Axis(fig[1, 2], xlabel = "field (T)", ylabel = "model - experiment (uB / Yb)",
           title = "Residuals. The disorder model has no plateau to explain away.")
hlines!(ax2, [0.0]; color = (:black, 0.4), linestyle = :dash)
for (k, r) in enumerate(results)
    lines!(ax2, Bs, r.obj.residual; color = cols[mod1(k, length(cols))], linewidth = 2.4,
           label = r.name)
end
for (k, r) in enumerate(results)
    lines!(ax2, Bs, r.model_nolin .- M_exp; color = cols[mod1(k, length(cols))],
           linewidth = 1.4, linestyle = :dot)
end
axislegend(ax2; position = :rt, labelsize = 10)
text!(ax2, Bs[1], 0.9 * maximum(vcat([r.obj.residual for r in results]...));
      text = "dotted: same set with the linear term forced to zero", fontsize = 10,
      align = (:left, :top), color = :grey30)

# ---- 3  dM/dH --------------------------------------------------------------------
ax3 = Axis(fig[2, 1], xlabel = "field (T)", ylabel = "dM/dH (uB / Yb / T)",
           title = "Differential susceptibility: a plateau would show as a dip toward zero")
lines!(ax3, Bs, deriv(Bs, M_exp); color = :black, linewidth = 2.0, label = "experiment")
scatter!(ax3, Bs, deriv(Bs, M_exp); color = :black, markersize = 6)
for (k, r) in enumerate(results)
    lines!(ax3, Bs, deriv(Bs, r.obj.model); color = cols[mod1(k, length(cols))],
           linewidth = 2.4, label = r.name)
end
axislegend(ax3; position = :rt, labelsize = 10)

# ---- 4  how much of the curve is the LINEAR term? --------------------------------
ax4 = Axis(fig[2, 2], xlabel = "field (T)", ylabel = "M (uB / Yb)",
           title = "Decomposition: is the fitted linear term physically plausible?")
scatter!(ax4, Bs, M_exp; color = :black, markersize = 7, label = "experiment")
for (k, r) in enumerate(results)
    r.name == "clean" && continue
    c = cols[mod1(k, length(cols))]
    spin = r.obj.A_M .* r.obj.raw
    lin = (r.obj.A_M * r.obj.chi_vv) .* Bs
    lines!(ax4, Bs, spin; color = c, linewidth = 2.2, label = "$(r.name): spin part")
    lines!(ax4, Bs, lin; color = c, linewidth = 2.0, linestyle = :dot,
           label = "$(r.name): linear part (chi_vv = $(round(r.obj.chi_vv; sigdigits=3)))")
end
# The crystal field gives chi_VV^zz = 0.0171 +- 0.0007 uB/T. Anything far above that is
# absorbing something other than Van Vleck, which is an open thread in its own right.
lines!(ax4, Bs, 0.0171 .* Bs; color = :grey40, linewidth = 2.0, linestyle = :dash,
       label = "crystal-field Van Vleck (0.0171 uB/T)")
axislegend(ax4; position = :lt, labelsize = 9)

Label(fig[4, 1:2],
      "Clean = same Hamiltonian with sigma_J = sigma_gzz = 0. It is a control isolating what " *
      "DISORDER does, not a reproduction of the published model, which uses XXZ with " *
      "Delta ~ 1.35 while this model is isotropic Heisenberg with all anisotropy in g.";
      fontsize = 11, color = :grey30)

# ---- 5  ABSOLUTE comparison, which profiling A_M out makes invisible -------------
ax5 = Axis(fig[3, 1], xlabel = "field (T)", ylabel = "M (uB / Yb)",
           title = "ABSOLUTE: no A_M, no fitted linear term (newly possible)")
scatter!(ax5, Bs, M_exp; color = :black, markersize = 7, label = "experiment (corrected MPMS3)")
for (k, r) in enumerate(results)
    c = cols[mod1(k, length(cols))]
    lines!(ax5, Bs, r.obj.raw; color = c, linewidth = 2.4, label = "$(r.name) raw spin")
end
# The crystal-field Van Vleck term is a genuine addition to the spin moment for an EFFECTIVE
# S = 1/2, so the honest absolute model is raw + chi_VV*B with chi_VV from the crystal field --
# NOT the fitted value, which the absolute data itself excludes.
lines!(ax5, Bs, results[1].obj.raw .+ 0.0171 .* Bs; color = :grey35, linewidth = 2.0,
       linestyle = :dash, label = "$(results[1].name) + CEF Van Vleck")
axislegend(ax5; position = :rb, labelsize = 8)

ax6 = Axis(fig[3, 2], xlabel = "field (T)", ylabel = "model / experiment",
           title = "Absolute ratio -- a SHAPE error, not a scale offset")
for (k, r) in enumerate(results)
    c = cols[mod1(k, length(cols))]
    ok = M_exp .> 0.05
    lines!(ax6, Bs[ok], (r.obj.raw .+ 0.0171 .* Bs)[ok] ./ M_exp[ok]; color = c, linewidth = 2.4,
           label = r.name)
end
hlines!(ax6, [1.0]; color = (:black, 0.55), linestyle = :dash)
axislegend(ax6; position = :rt, labelsize = 9)
text!(ax6, 0.05, 0.06; space = :relative, fontsize = 9.5, color = :grey25,
      text = "A ratio that DECLINES with field is a shape error: the model polarises too
" *
             "easily at low field and the two converge only as saturation is forced.
" *
             "A pure normalisation error would be a flat offset. This is what profiling
" *
             "A_M out was hiding, and it is the defect the spurious linear term patched.")

out = joinpath(FDIR, "mvh_parameter_sets.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
println("wrote " * joinpath(TDIR, "mvh_parameter_sets.csv"))
