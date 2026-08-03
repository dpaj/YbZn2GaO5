#!/usr/bin/env julia
# Overplot the six 1D neutron cuts against experiment for SEVERAL NAMED PARAMETER SETS.
#
# Distinct from plot_neutron_vs_exp.jl, which varies sigma_J alone at fixed everything
# else. This one compares whole parameter sets -- by-eye against a fitted set -- which is
# the comparison that shows whether a fit actually improved the LINESHAPE or merely the
# chi2. One global intensity scale is profiled out per set by weighted least squares, so
# only lineshape is being compared; Sunny's prefactor is not comparable to the analytical
# model's and the experimental normalisation is itself suspect.
#
#   julia -t auto --project=. scripts/plot_neutron_parameter_sets.jl
#
# Parameter sets come from [[sets]] in the config. Cost is one full neutron evaluation per
# set, so keep the list short.

using Printf, Statistics, LinearAlgebra, CairoMakie, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const LOADED = SV.sv_load_diagnostic_controls(REPO_ROOT,
    "configs/neutron_parameter_sets_controls.toml";
    env_var="YZGO_NEUTRON_SETS_CONTROLS")
const CFG = LOADED.diag
const RUN = get(CFG, "run", Dict{String,Any}())
const controls = LOADED.controls

params0, applied = SV.sv_apply_param_overrides(
    SV.sv_load_params(REPO_ROOT, controls).params, RUN)

const CUTS = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
const NREAL = Int(get(RUN, "n_realizations", 4))
const MAXITERS = Int(get(RUN, "minimize_maxiters", 1000))
const RELAX = Int(get(RUN, "relax_attempts", 1))
const ELASTIC_CUT = Float64(get(RUN, "elastic_cutoff_meV", 0.35))
const FDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "figure_subdir",
    "results/figures/sunny_validation/neutron_parameter_sets"))
const TDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "table_subdir",
    "results/feature_tables/sunny_validation/neutron_parameter_sets"))
mkpath(FDIR); mkpath(TDIR)

"Pretty qtag: 0p33_0p33_0 -> (0.33, 0.33, 0)"
function qlabel(t::AbstractString)
    p = split(t, "_")
    length(p) == 3 || return t
    return "(" * join(replace.(p, "p" => "."), ", ") * ")"
end

const PLOT_ONLY = Bool(get(RUN, "plot_only", false)) ||
                 !isempty(get(ENV, "YZGO_PLOT_ONLY", ""))
sets = get(CFG, "sets", Any[])
(PLOT_ONLY || !isempty(sets)) || error("No [[sets]] in the config; nothing to compare.")

@printf("threads %d, realizations %d, cuts %d, q/cut %d, quadrature %s\n",
        Threads.nthreads(), NREAL, length(CUTS),
        length(SV.sv_kpm_1d_q_sampler(CUTS[1], controls).qs),
        string(get(controls["kpm"]["q_averaging"], "resolution_quadrature", "gauss_hermite")))
isempty(applied) || @printf("base overrides: %s\n", join(applied, ", "))

results = NamedTuple[]
for (i, sd) in (PLOT_ONLY ? () : enumerate(sets))
    name = String(get(sd, "name", "set$i"))
    skip = ("name", "n_realizations")
    over = NamedTuple(Symbol(k) => Float64(v) for (k, v) in sd if !(k in skip))
    p = merge(params0, over)
    # A clean set (sigma_J = sigma_gzz = 0) draws the SAME system for every realization, so
    # averaging over more than one is exactly redundant. Let a set override the count.
    nr = Int(get(sd, "n_realizations", NREAL))
    (p.sigma_J == 0 && p.sigma_gzz == 0 && nr > 1) &&
        @warn "set $name has no disorder; its realizations are identical, so n_realizations > 1 is wasted" nr
    t = time()
    o = SV.sv_neutron_objective(p, controls, CUTS; realizations=0:(nr - 1),
        threaded=true, maxiters=MAXITERS, relax_attempts=RELAX, on_failure=:record)
    @printf("\n%-12s J1=%.3f sigma_J=%.2f gzz=%.2f sigma_gzz=%.2f -> chi2_red = %.4g  (%.0f s)\n",
            name, p.J1_meV, p.sigma_J, p.gzz, p.sigma_gzz, o.chi2_red, time() - t)
    for pc in o.per_cut
        @printf("    %-14s %5.1f T  chi2_red = %-10.4g rms = %.4g\n",
                pc.qtag, pc.field_T, pc.chi2_red, pc.rms)
    end
    push!(results, (; name, params=p, obj=o))
end

# ---------------------------------------------------------------- persist
# Written BEFORE any plotting, and the figure is then built by reading these back. A
# plotting bug therefore costs a rerun of the plot, never a rerun of the KPM -- which it
# did cost once, to a @sprintf format-string error after 950 s of compute.
const CHI2_CSV = joinpath(TDIR, "per_cut_chi2.csv")
const CURVE_CSV = joinpath(TDIR, "model_curves.csv")

if !PLOT_ONLY
    open(CHI2_CSV, "w") do io
        println(io, "set,J1_meV,sigma_J,gzz,sigma_gzz,chi2_red_total,scale,field_T,qtag,chi2_red,rms")
        for r in results, pc in r.obj.per_cut
            @printf(io, "%s,%.4f,%.4f,%.4f,%.4f,%.6g,%.6g,%.1f,%s,%.6g,%.6g\n",
                    r.name, r.params.J1_meV, r.params.sigma_J, r.params.gzz,
                    r.params.sigma_gzz, r.obj.chi2_red, r.obj.scale, pc.field_T, pc.qtag,
                    pc.chi2_red, pc.rms)
        end
    end
    open(CURVE_CSV, "w") do io
        println(io, "set,qtag,field_T,energy_meV,I_exp,Ierr_exp,I_model_scaled,cut_chi2_red")
        for r in results
            for (i, cut) in enumerate(CUTS)
                m = r.obj.scale .* r.obj.curves[i]
                c2 = r.obj.per_cut[i].chi2_red
                for j in eachindex(cut.energy_meV)
                    @printf(io, "%s,%s,%.1f,%.6g,%.6g,%.6g,%.6g,%.6g\n",
                            r.name, cut.qtag, cut.field_T, cut.energy_meV[j],
                            cut.intensity[j], cut.error[j], m[j], c2)
                end
            end
        end
    end
    println("\nwrote $CHI2_CSV")
    println("wrote $CURVE_CSV")
end

# ---------------------------------------------------------------- figure (from CSV)
isfile(CURVE_CSV) || error("No $CURVE_CSV; run once without plot_only to generate it.")
hdr = String.(split(strip(readlines(CURVE_CSV)[1]), ','))
crows = Dict{String,String}[]
for l in readlines(CURVE_CSV)[2:end]
    isempty(strip(l)) && continue
    f = String.(split(strip(l), ','))
    length(f) == length(hdr) && push!(crows, Dict(zip(hdr, f)))
end
cnum(r, k) = something(tryparse(Float64, get(r, k, "")), NaN)

setnames = unique(String[r["set"] for r in crows])
qtags = unique(String[r["qtag"] for r in crows])
fields = sort(unique(Float64[cnum(r, "field_T") for r in crows]))

fig = Figure(size = (520 * length(qtags), 430 * length(fields)))
Label(fig[0, 1:length(qtags)],
      "YbZn2GaO5 -- 1D neutron cuts vs model, $(NREAL) realizations, " *
      "one global intensity scale profiled out per parameter set";
      fontsize = 17, font = :bold)
cols = [:crimson, :dodgerblue, :seagreen, :darkorange, :purple]

for (row, B) in enumerate(fields), (col, qt) in enumerate(qtags)
    sel(nm) = [r for r in crows if r["qtag"] == qt && cnum(r, "field_T") ≈ B &&
                                   r["set"] == nm]
    first_rows = sel(setnames[1])
    isempty(first_rows) && continue
    ax = Axis(fig[row, col], xlabel = "energy transfer (meV)",
              ylabel = "intensity (arb.)",
              title = "$(qlabel(qt))   $(round(Int, B)) T")
    E = [cnum(r, "energy_meV") for r in first_rows]
    y = [cnum(r, "I_exp") for r in first_rows]
    e = [cnum(r, "Ierr_exp") for r in first_rows]
    errorbars!(ax, E, y, e; color = (:black, 0.45), whiskerwidth = 4)
    scatter!(ax, E, y; color = :black, markersize = 7, label = "experiment")
    # A clean model's magnon is limited only by the instrument and the KPM kernel, so it is
    # far taller and narrower than the data. Letting it set the y range would flatten every
    # other curve, so it is excluded from the scale and allowed to clip.
    noscale = String.(get(RUN, "ylim_exclude_sets", ["clean"]))
    vals = Float64[]
    append!(vals, filter(isfinite, y[E .>= ELASTIC_CUT]))
    for (k, nm) in enumerate(setnames)
        rs = sel(nm)
        isempty(rs) && continue
        Em = [cnum(r, "energy_meV") for r in rs]
        m = [cnum(r, "I_model_scaled") for r in rs]
        c2 = cnum(rs[1], "cut_chi2_red")
        lines!(ax, Em, m; color = cols[mod1(k, length(cols))], linewidth = 2.6,
               label = "$nm  (chi2=$(round(Int, c2)))")
        nm in noscale || append!(vals, filter(isfinite, m[Em .>= ELASTIC_CUT]))
    end
    # The elastic line and any quasi-elastic divergence are orders of magnitude above the
    # magnon signal and would flatten every panel, so scale y from the inelastic region.
    if !isempty(vals)
        hi = maximum(vals)
        ylims!(ax, -0.06 * hi, 1.20 * hi)
    end
    xlims!(ax, 0.0, maximum(E))
    vspan!(ax, 0.5, 3.0; color = (:grey, 0.10))
    row == 1 && col == 1 && axislegend(ax; position = :rt, labelsize = 10)
end
Label(fig[length(fields) + 1, 1:length(qtags)],
      "Shaded band is the [0.5, 3.0] meV fit window. Outside it the model is polluted by " *
      "Goldstone-like divergence as hw -> 0, and the data by elastic, incoherent and " *
      "possibly Bragg scattering, so neither side is meaningful there.";
      fontsize = 11, color = :grey30)

out = joinpath(FDIR, "neutron_parameter_sets.png")
save(out, fig; px_per_unit = 2)
println("wrote $out")
