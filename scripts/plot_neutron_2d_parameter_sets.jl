#!/usr/bin/env julia
# Side-by-side 2D data-vs-model comparison along the CNCS path, for one or more NAMED
# parameter sets. NOT a fit -- a forward comparison at parameters chosen elsewhere.
#
#   julia -t auto --project=. scripts/plot_neutron_2d_parameter_sets.jl
#   YZGO_PLOT_ONLY=1 YZGO_PLOT_SETS="fitted" julia --project=. scripts/plot_neutron_2d_parameter_sets.jl
#
# This replaces the legacy sv_run_kpm_2d_data_model_comparison route for model computation.
# That path had drifted from the 1D one in two ways that matter:
#
#   * it started from `randomize_spins!` rather than a field-polarized state, so in a
#     disordered system at partial saturation it relaxed into an arbitrary local minimum --
#     the 2D map was not reproducible and was not the same state the 1D cuts came from;
#   * it constructed SpinWaveTheoryKPM without passing `regularization`, so
#     [kpm].regularization was silently ignored and Sunny's 1e-8 default used regardless.
#
# Going through sv_neutron_2d_curves fixes both and additionally gets q-threading, the
# pooled KPM operators, relax_attempts = 1 and realization averaging, none of which the
# legacy path had.
#
# One intensity scale is fitted per parameter set by nonnegative least squares over both
# fields, and the SAME colour range (taken from the data) is used for every panel. Panels on
# different colour scales would make a side-by-side comparison meaningless.

using Printf, Statistics, LinearAlgebra, CairoMakie, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const LOADED = SV.sv_load_diagnostic_controls(REPO_ROOT,
    "configs/neutron_2d_parameter_sets_controls.toml"; env_var="YZGO_2D_SETS_CONTROLS")
const CFG = LOADED.diag
const RUN = get(CFG, "run", Dict{String,Any}())
const controls = LOADED.controls

params0, applied = SV.sv_apply_param_overrides(
    SV.sv_load_params(REPO_ROOT, controls).params, RUN)

const K2 = SV.sv_kpm_2d_controls(controls)
const LEG = Int(get(K2, "leg", 1))
const FIELDS = Float64.(get(K2, "fields_T", [9.0, 14.0]))
const SCANS = SV.sv_load_2d_scans_for_kpm(REPO_ROOT, controls; fields_T=FIELDS, leg=LEG)
isempty(SCANS) && error("No experimental 2D scans loaded.")

const NREAL = Int(get(RUN, "n_realizations", 2))
const MAXITERS = Int(get(RUN, "minimize_maxiters", 1000))
const RELAX = Int(get(RUN, "relax_attempts", 1))
const FDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "figure_subdir",
    "results/figures/sunny_validation/neutron_2d_parameter_sets"))
const TDIR = SV.sv_repo_path(REPO_ROOT, get(controls["paths"], "table_subdir",
    "results/feature_tables/sunny_validation/neutron_2d_parameter_sets"))
mkpath(FDIR); mkpath(TDIR)
const CURVE_CSV = joinpath(TDIR, "model_maps.csv")
const PLOT_ONLY = Bool(get(RUN, "plot_only", false)) ||
                  !isempty(get(ENV, "YZGO_PLOT_ONLY", ""))

sets = get(CFG, "sets", Any[])
(PLOT_ONLY || !isempty(sets)) || error("No [[sets]] in the config.")

if !PLOT_ONLY
    nx = length(first(values(SCANS)).x)
    sampler = SV.sv_kpm_2d_q_sampler_from_scan(first(values(SCANS)), K2; leg=LEG)
    @printf("threads %d, fields %s, leg %d\n", Threads.nthreads(), string(FIELDS), LEG)
    @printf("path points %d, q evaluated per field %d (%.1f per pixel), energies %d -> %d\n",
            nx, length(sampler.qs_flat), length(sampler.qs_flat) / nx,
            Int(controls["kpm"]["n_energy"]), length(first(values(SCANS)).e))
    @printf("quadrature %s, tol %s, regularization %s\n",
            string(get(get(controls["kpm_2d"], "q_averaging", Dict()), "resolution_quadrature", "gauss_hermite")),
            string(controls["kpm"]["tol"]), string(get(controls["kpm"], "regularization", "?")))
    isempty(applied) || @printf("base overrides: %s\n", join(applied, ", "))
end

results = NamedTuple[]
for (i, sd) in (PLOT_ONLY ? () : enumerate(sets))
    name = String(get(sd, "name", "set$i"))
    skip = ("name", "n_realizations")
    over = NamedTuple(Symbol(k) => Float64(v) for (k, v) in sd if !(k in skip))
    p = merge(params0, over)
    nr = Int(get(sd, "n_realizations", NREAL))
    (p.sigma_J == 0 && p.sigma_gzz == 0 && nr > 1) &&
        @warn "set $name has no disorder, so its realizations are identical" nr
    t = time()
    r = SV.sv_neutron_2d_curves(p, controls, SCANS, K2; realizations=0:(nr - 1), leg=LEG,
            threaded=true, maxiters=MAXITERS, relax_attempts=RELAX, verbose=true)
    # Deposit onto each scan's own energy grid, then fit ONE scale across both fields.
    on_grid = Dict{Float64,Matrix{Float64}}()
    for B in keys(r.curves)
        on_grid[B] = SV.sv_model_to_scan_energy_grid(r.energy_meV, r.curves[B], SCANS[B];
                         controls=controls, section="kpm_2d")
    end
    scale = SV.sv_least_squares_scale_2d(SCANS, on_grid;
                emin=Float64(get(K2, "scale_energy_min_meV", 0.25)),
                emax=Float64(get(K2, "scale_energy_max_meV", 3.2)), fallback=1.0)
    @printf("%-16s J1=%.3f sigma_J=%.2f gzz=%.2f sigma_gzz=%.2f  scale=%.5g  reg=%s  (%.0f s, KPM %.0f s)\n",
            name, p.J1_meV, p.sigma_J, p.gzz, p.sigma_gzz, scale,
            string(r.regularization_values), time() - t, r.kpm_seconds)
    push!(results, (; name, params=p, n_real=nr, on_grid, scale))
    flush(stdout)
end

# ---------------------------------------------------------------- persist
# Written before plotting so a plotting bug never costs the KPM again, and the figure is
# built by reading this back, so the fresh and replay paths are the same code.
if !PLOT_ONLY
    open(CURVE_CSV, "w") do io
        println(io, "set,J1_meV,sigma_J,gzz,sigma_gzz,n_realizations,scale,field_T," *
                    "path_coordinate,energy_meV,I_data,I_model_scaled")
        for r in results, B in sort(collect(keys(r.on_grid)))
            scan = SCANS[B]; m = r.on_grid[B]
            for ix in eachindex(scan.x), ie in eachindex(scan.e)
                @printf(io, "%s,%.4f,%.4f,%.4f,%.4f,%d,%.6g,%.1f,%.6g,%.6g,%.6g,%.6g\n",
                        r.name, r.params.J1_meV, r.params.sigma_J, r.params.gzz,
                        r.params.sigma_gzz, r.n_real, r.scale, B, scan.x[ix], scan.e[ie],
                        scan.z[ix, ie], r.scale * m[ix, ie])
            end
        end
    end
    println("wrote $CURVE_CSV")
end

# ---------------------------------------------------------------- figure (from CSV)
isfile(CURVE_CSV) || error("No $CURVE_CSV; run once without plot_only to generate it.")
hdr = String.(split(strip(readlines(CURVE_CSV)[1]), ','))
rows = Dict{String,String}[]
for l in readlines(CURVE_CSV)[2:end]
    isempty(strip(l)) && continue
    f = String.(split(strip(l), ','))
    length(f) == length(hdr) && push!(rows, Dict(zip(hdr, f)))
end
gnum(r, k) = something(tryparse(Float64, get(r, k, "")), NaN)

allsets = unique(String[r["set"] for r in rows])
wanted = let e = get(ENV, "YZGO_PLOT_SETS", "")
    isempty(e) ? String.(get(RUN, "plot_sets", String[])) : String.(strip.(split(e, ',')))
end
setnames = isempty(wanted) ? allsets : [s for s in allsets if s in wanted]
isempty(setnames) && error("plot_sets $(wanted) matched none of $(allsets)")
fieldvals = sort(unique(Float64[gnum(r, "field_T") for r in rows]))

for nm in setnames
    r = first(f for f in rows if f["set"] == nm)
    @printf("  %-16s J1=%s sigma_J=%s gzz=%s sigma_gzz=%s  (%s realizations, scale %s)\n",
            nm, r["J1_meV"], r["sigma_J"], r["gzz"], r["sigma_gzz"],
            r["n_realizations"], r["scale"])
end

"Reassemble an (x, e) matrix for one set and field from the long CSV."
function grid_of(setname, B; key="I_model_scaled")
    sel = [r for r in rows if r["set"] == setname && gnum(r, "field_T") ≈ B]
    xs = sort(unique(gnum.(sel, "path_coordinate")))
    es = sort(unique(gnum.(sel, "energy_meV")))
    xi = Dict(v => i for (i, v) in enumerate(xs)); ei = Dict(v => i for (i, v) in enumerate(es))
    Z = fill(NaN, length(xs), length(es))
    for r in sel
        Z[xi[gnum(r, "path_coordinate")], ei[gnum(r, "energy_meV")]] = gnum(r, key)
    end
    return (xs, es, Z)
end

ylim = Float64.(get(K2, "energy_ylim_meV", [0.20, 3.20]))
guides = Float64.(get(K2, "guide_xs", Float64[]))
tickpos = Float64.(get(K2, "xtick_positions", Float64[]))
ticklab = String.(get(K2, "xtick_labels", String[]))

# One colour range for EVERY panel, taken from the data. Different scales per panel would
# make the comparison meaningless, which is the whole point of the figure.
dat_ranges = Float64[]
for B in fieldvals
    _, es, Z = grid_of(setnames[1], B; key="I_data")
    sel = (es .>= Float64(get(K2, "data_color_energy_min_meV", 0.25)))
    append!(dat_ranges, filter(isfinite, vec(Z[:, sel])))
end
hi = isempty(dat_ranges) ? 1.0 :
     quantile(dat_ranges, Float64(get(K2, "data_clip_high_quantile", 0.95)))
crange = (0.0, hi)

nrow = 1 + length(setnames)
fig = Figure(size = (560 * length(fieldvals) + 120, 330 * nrow + 90))
Label(fig[0, 1:length(fieldvals)],
      "YbZn2GaO5 -- 2D CNCS path, data vs model (forward comparison, not a fit). " *
      "One fitted intensity scale per set; identical colour range on every panel.";
      fontsize = 16, font = :bold)

for (col, B) in enumerate(fieldvals)
    xs, es, D = grid_of(setnames[1], B; key="I_data")
    ax = Axis(fig[1, col], ylabel = col == 1 ? "energy (meV)" : "",
              title = @sprintf("experiment   %.0f T", B))
    heatmap!(ax, xs, es, D; colorrange = crange,
             colormap = Symbol(get(K2, "colormap", "viridis")))
    ylims!(ax, ylim[1], ylim[2])
    isempty(guides) || vlines!(ax, guides; color = (:white, 0.35), linewidth = 1)
    isempty(tickpos) || (ax.xticks = (tickpos, ticklab))
    for (irow, nm) in enumerate(setnames)
        xs2, es2, M = grid_of(nm, B)
        ax2 = Axis(fig[1 + irow, col], ylabel = col == 1 ? "energy (meV)" : "",
                   xlabel = irow == length(setnames) ? "path coordinate" : "",
                   title = @sprintf("%s   %.0f T", nm, B))
        hm = heatmap!(ax2, xs2, es2, M; colorrange = crange,
                      colormap = Symbol(get(K2, "colormap", "viridis")))
        ylims!(ax2, ylim[1], ylim[2])
        isempty(guides) || vlines!(ax2, guides; color = (:white, 0.35), linewidth = 1)
        isempty(tickpos) || (ax2.xticks = (tickpos, ticklab))
        (col == length(fieldvals) && irow == length(setnames)) &&
            Colorbar(fig[2:nrow, length(fieldvals) + 1], hm, label = "intensity (arb.)")
    end
end

Label(fig[nrow + 1, 1:length(fieldvals)],
      "Model computed from a FIELD-POLARIZED relaxed state with [kpm].regularization " *
      "honoured; the legacy 2D route used a random start and silently ignored " *
      "regularization. Energy resolution deposited onto the experimental grid.";
      fontsize = 10, color = :grey30)

suffix = length(setnames) == length(allsets) ? "" : "_" * join(setnames, "_")
out = joinpath(FDIR, "neutron_2d_parameter_sets$(suffix).png")
save(out, fig; px_per_unit = 2)
println("wrote $out")
