#!/usr/bin/env julia
# Overplot the Ei = 3.32 meV and Ei = 4.65 meV 1D energy scans. DATA ONLY -- no model, no
# Sunny, no KPM, so it runs in about a minute.
#
#   julia --project=. scripts/plot_neutron_1d_two_incident_energies.jl
#
# WHY THIS IS THE RIGHT BACKGROUND DIAGNOSTIC
#
# The sample signal at a given (Q, E) is a property of the sample and must be the same at
# both incident energies. Instrumental background is NOT: it depends on the spectrometer
# configuration, so the same (Q, E) picks up a different background contribution at
# different Ei. So a feature that survives at the same energy with the same relative
# strength in both datasets is sample signal, and one that changes is not.
#
# That matters because the sharp, nearly flat but slightly curved feature near 2 meV is a
# MAGNET background. It is what the background construction exists to remove, it is present
# at both 9 and 14 T, and at 14 T it overlaps the real mode position almost exactly -- so at
# 14 T the two cannot be separated by eye within a single dataset.
#
# `(0,1,0)` is NOT kinematically accessible at Ei = 3.32, which is unlucky: it is also the
# cut with no structured-residual correction, so it is the one cut that can neither be
# cross-checked here nor corrected there.
#
# Everything plotted is RAW, before any background subtraction, which is the point.
#
# The relative normalisation between the two Ei is FITTED (different flux, different
# detector coverage), so only shapes and positions are comparable, never absolute heights.
# Note also that Ei = 3.32 has BETTER energy resolution, so a genuinely sharp feature should
# look sharper there -- that is a discriminator, not a nuisance.

using Printf, Statistics, LinearAlgebra, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const controls = SV.sv_load_controls(REPO_ROOT)
const DIR = SV.sv_repo_path(REPO_ROOT, controls["paths"]["neutron_1d_dir"])
const FDIR = SV.sv_repo_path(REPO_ROOT,
    "results/figures/sunny_validation/two_incident_energies")
mkpath(FDIR)

# Load every 1D scan present, keyed by (Ei, T, field, qtag).
scans = Dict{Tuple{Float64,Float64,Float64,String},Any}()
for f in sort(readdir(DIR))
    endswith(lowercase(f), ".dat") || continue
    path = joinpath(DIR, f)
    meta = try
        SV.sv_parse_neutron_1d_filename(path)
    catch
        continue
    end
    scans[(meta.Ei_meV, meta.temperature_K, meta.field_T, meta.qtag)] =
        SV.sv_load_neutron_raw_scan_1d(path, controls)
end
@printf("loaded %d raw 1D scans\n", length(scans))

const EI_HI = 4.65
const EI_LO = 3.32
const TBASE = 0.07

qtags = sort(unique(String[k[4] for k in keys(scans)]))
fields = sort(unique(Float64[k[3] for k in keys(scans)]))
# Only qtags measured at BOTH incident energies are informative here.
both = [q for q in qtags if any(k -> k[1] == EI_LO && k[4] == q, keys(scans)) &&
                            any(k -> k[1] == EI_HI && k[4] == q, keys(scans))]
missing_lo = [q for q in qtags if !(q in both)]
@printf("qtags at both Ei: %s\n", join(both, ", "))
isempty(missing_lo) || @printf("qtags only at Ei = %.2f: %s  (not accessible at %.2f)\n",
                               EI_HI, join(missing_lo, ", "), EI_LO)

"Pretty qtag: 0p33_0p33_0 -> (0.33, 0.33, 0)"
qlabel(t) = (p = split(t, "_"); length(p) == 3 ?
             "(" * join(replace.(p, "p" => "."), ", ") * ")" : t)

"""
Single nonnegative scale putting the low-Ei scan on the high-Ei intensity scale, fitted over
the overlapping energy range above `emin`. Different flux and detector coverage mean the two
are not absolutely comparable; this makes SHAPES comparable and nothing more.
"""
function relative_scale(lo, hi; emin=0.35)
    yl = SV.sv_interp1(lo.energy_meV, lo.intensity, hi.energy_meV)
    ok = isfinite.(yl) .& isfinite.(hi.intensity) .& (hi.energy_meV .>= emin)
    count(ok) < 5 && return 1.0
    num = dot(yl[ok], hi.intensity[ok]); den = dot(yl[ok], yl[ok])
    return (isfinite(den) && den > 0) ? max(0.0, num / den) : 1.0
end

fig = Figure(size = (560 * length(both) + 60, 340 * length(fields) + 110))
Label(fig[0, 1:length(both)],
      "YbZn2GaO5 -- RAW 1D scans at two incident energies. Sample signal must appear at the " *
      "same energy in both; instrumental background need not.";
      fontsize = 16, font = :bold)

anchor_lo = Float64.(get(controls["kpm"], "min_bg_low_window_meV", [0.0, 0.75]))
anchor_hi = Float64(get(controls["kpm"], "min_bg_high_threshold_meV", 2.5))

for (row, B) in enumerate(fields), (col, q) in enumerate(both)
    hi = get(scans, (EI_HI, TBASE, B, q), nothing)
    lo = get(scans, (EI_LO, TBASE, B, q), nothing)
    (hi === nothing && lo === nothing) && continue
    ax = Axis(fig[row, col], xlabel = row == length(fields) ? "energy transfer (meV)" : "",
              ylabel = col == 1 ? "raw intensity (arb.)" : "",
              title = @sprintf("%s   %.0f T", qlabel(q), B))

    # Shade the background construction's footing: it anchors on min-over-fields inside the
    # two windows and PCHIP-INTERPOLATES everything between, which is ~70% of the [0.5, 3.0]
    # fit window and is exactly where the magnet feature lives.
    vspan!(ax, anchor_lo[1], anchor_lo[2]; color = (:seagreen, 0.13))
    vspan!(ax, anchor_hi, 4.2; color = (:seagreen, 0.13))
    vspan!(ax, anchor_lo[2], anchor_hi; color = (:indianred, 0.09))

    s = (hi !== nothing && lo !== nothing) ? relative_scale(lo, hi) : 1.0
    if hi !== nothing
        errorbars!(ax, hi.energy_meV, hi.intensity, hi.error; color = (:black, 0.3),
                   whiskerwidth = 3)
        scatter!(ax, hi.energy_meV, hi.intensity; color = :black, markersize = 6,
                 label = @sprintf("Ei = %.2f meV", EI_HI))
        lines!(ax, hi.energy_meV, hi.intensity; color = (:black, 0.5), linewidth = 1)
    end
    if lo !== nothing
        errorbars!(ax, lo.energy_meV, s .* lo.intensity, s .* lo.error;
                   color = (:dodgerblue, 0.3), whiskerwidth = 3)
        scatter!(ax, lo.energy_meV, s .* lo.intensity; color = :dodgerblue, markersize = 6,
                 label = @sprintf("Ei = %.2f meV (x%.3g)", EI_LO, s))
        lines!(ax, lo.energy_meV, s .* lo.intensity; color = (:dodgerblue, 0.5), linewidth = 1)
    end
    # The 20 K scan exists only at Ei = 3.32 and 0 T. At 20 K the magnetic signal is largely
    # thermally destroyed, so what remains is mostly background -- an independent probe of it.
    hot = get(scans, (EI_LO, 20.0, B, q), nothing)
    if hot !== nothing
        sh = hi !== nothing ? relative_scale(hot, hi) : 1.0
        lines!(ax, hot.energy_meV, sh .* hot.intensity; color = :darkorange, linewidth = 2.2,
               linestyle = :dash, label = @sprintf("20 K, Ei = %.2f (x%.3g)", EI_LO, sh))
    end

    vals = Float64[]
    hi === nothing || append!(vals, filter(isfinite, hi.intensity[hi.energy_meV .>= 0.35]))
    lo === nothing || append!(vals, filter(isfinite, s .* lo.intensity[lo.energy_meV .>= 0.35]))
    isempty(vals) || ylims!(ax, -0.08 * maximum(vals), 1.25 * maximum(vals))
    xlims!(ax, 0.0, 4.0)
    row == 1 && col == 1 && axislegend(ax; position = :rt, labelsize = 10)
end

Label(fig[length(fields) + 1, 1:length(both)],
      "Green = the two windows where the background is ANCHORED on min-over-fields data. " *
      "Red = the gap it PCHIP-interpolates across, about 70% of the [0.5, 3.0] meV fit " *
      "window, and where the magnet feature sits. Relative normalisation between the two Ei " *
      "is fitted, so compare positions and shapes only, never absolute heights. Ei = 3.32 " *
      "has the better energy resolution, so a genuinely sharp feature should look sharper " *
      "there.";
      fontsize = 10, color = :grey30)

out = joinpath(FDIR, "neutron_1d_two_incident_energies.png")
save(out, fig; px_per_unit = 2)
println("wrote $out")
