#!/usr/bin/env julia
# Overplot the Ei = 3.32 meV and Ei = 4.65 meV 1D energy scans. DATA ONLY -- no model, no
# Sunny, no KPM, so it runs in about a minute.
#
#   julia --project=. scripts/plot_neutron_1d_two_incident_energies.jl
#   YZGO_EI_OVERLAY=full julia --project=. scripts/plot_neutron_1d_two_incident_energies.jl
#
# TWO MODES. The default ("raw") is the clean two-dataset comparison and is deliberately
# left alone -- it is the figure that isolates the background question with nothing else in
# the way. "full" adds the background that was actually subtracted, the corrected data, the
# current best-fit model, and a raw estimate of the Ei = 4.65 background from the difference
# between the two datasets. Both write separate files, so neither overwrites the other.
#
# "full" also includes `(0,1,0)`, which has no Ei = 3.32 counterpart, because its 0 T scan is
# still a background probe. NOTE the 0 T data are NOT signal-free -- there is still a magnon,
# just at lower energy, plus a diffuse zero-field continuum. What makes 0 T useful is that it
# is expected signal-free in a LIMITED HIGH-ENERGY region, so the absence of a 2.08 meV spike
# there is informative about background at that Q even though the scan as a whole is not a
# background measurement.
#
# This is the same logic the background construction itself uses, and it is worth stating
# because it explains the interpolation gap rather than merely describing it. The magnet
# background sits at a fixed energy while the SIGNAL MOVES WITH FIELD, so min-over-fields
# picks, at each energy, whichever field has its signal furthest away. The 0.75-2.5 meV gap
# is precisely the range the signal SWEEPS THROUGH between 0 and 14 T, so no field leaves it
# clean -- hence interpolation. By the same token the LOW-energy side of the background is
# best estimated at 14 T, where the signal has moved up to ~2.1 meV and vacated it.
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
const MODE = lowercase(get(ENV, "YZGO_EI_OVERLAY", "raw"))
MODE in ("raw", "full") || error("YZGO_EI_OVERLAY must be raw or full, got $MODE")

# In "full" mode, pull in the background-corrected cuts (which carry raw_intensity,
# background and corrected intensity on ONE consistent scale, since intensity = raw - bg)
# and the saved best-fit model curves. Nothing is recomputed: the model comes from the CSV
# that scripts/plot_neutron_parameter_sets.jl already wrote, so this stays a data-only script.
corrected = Dict{Tuple{Float64,String},Any}()
model = Dict{Tuple{Float64,String},Tuple{Vector{Float64},Vector{Float64}}}()
if MODE == "full"
    for c in SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
        corrected[(c.field_T, c.qtag)] = c
    end
    mp = SV.sv_repo_path(REPO_ROOT,
        "results/feature_tables/sunny_validation/neutron_parameter_sets/model_curves.csv")
    if isfile(mp)
        lines = filter(!isempty, strip.(readlines(mp)))
        hdr = String.(split(lines[1], ','))
        want = get(ENV, "YZGO_MODEL_SET", "fitted")
        acc = Dict{Tuple{Float64,String},Vector{Tuple{Float64,Float64}}}()
        for l in lines[2:end]
            f = String.(split(l, ','))
            length(f) == length(hdr) || continue
            d = Dict(zip(hdr, f))
            get(d, "set", "") == want || continue
            B = something(tryparse(Float64, get(d, "field_T", "")), NaN)
            E = something(tryparse(Float64, get(d, "energy_meV", "")), NaN)
            I = something(tryparse(Float64, get(d, "I_model_scaled", "")), NaN)
            (isfinite(B) && isfinite(E)) || continue
            push!(get!(acc, (B, get(d, "qtag", "")), Tuple{Float64,Float64}[]), (E, I))
        end
        for (k, v) in acc
            sort!(v; by = first)
            model[k] = (first.(v), last.(v))
        end
        @printf("loaded model set '%s' for %d (field, qtag) pairs\n", want, length(model))
    else
        @warn "No saved model curves; run scripts/plot_neutron_parameter_sets.jl first" path=mp
    end
end

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
function relative_scale(lo, hi; emin=0.35, wins=nothing)
    yl = SV.sv_interp1(lo.energy_meV, lo.intensity, hi.energy_meV)
    inwin = wins === nothing ? (hi.energy_meV .>= emin) :
            reduce((a, b) -> a .| b,
                   [(hi.energy_meV .>= w[1]) .& (hi.energy_meV .<= w[2]) for w in wins])
    ok = isfinite.(yl) .& isfinite.(hi.intensity) .& inwin
    count(ok) < 5 && return 1.0
    num = dot(yl[ok], hi.intensity[ok]); den = dot(yl[ok], yl[ok])
    return (isfinite(den) && den > 0) ? max(0.0, num / den) : 1.0
end

# "raw" compares only qtags present at BOTH energies, which is the whole point of that
# figure. "full" shows every qtag, because (0,1,0) still has a 0 T background probe.
cols_q = MODE == "full" ? qtags : both

fig = Figure(size = (560 * length(cols_q) + 60, 340 * length(fields) + 130))
Label(fig[0, 1:length(cols_q)],
      MODE == "full" ?
      "YbZn2GaO5 -- 1D scans: raw at both incident energies, the subtracted background, the corrected data, and the best-fit model" :
      "YbZn2GaO5 -- RAW 1D scans at two incident energies. Sample signal must appear at the same energy in both; instrumental background need not.";
      fontsize = 16, font = :bold)

anchor_lo = Float64.(get(controls["kpm"], "min_bg_low_window_meV", [0.0, 0.75]))
anchor_hi = Float64(get(controls["kpm"], "min_bg_high_threshold_meV", 2.5))

for (row, B) in enumerate(fields), (col, q) in enumerate(cols_q)
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

    # In "full" mode anchor the relative normalisation on the windows where the background
    # construction has real data, and NOT across the interpolated gap. Then the difference
    # between the two datasets inside the gap is readable as background rather than being
    # partly absorbed into the scale.
    s = if hi !== nothing && lo !== nothing
        MODE == "full" ?
            relative_scale(lo, hi; wins=[(0.35, anchor_lo[2]), (anchor_hi, 2.9)]) :
            relative_scale(lo, hi)
    else
        1.0
    end
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

    if MODE == "full"
        c = get(corrected, (B, q), nothing)
        if c !== nothing
            # intensity = raw_intensity - background, all on one scale.
            lines!(ax, c.energy_meV, c.background; color = :darkorange, linewidth = 2.2,
                   linestyle = :dashdot, label = "background subtracted")
            errorbars!(ax, c.energy_meV, c.intensity, c.error; color = (:seagreen, 0.35),
                       whiskerwidth = 3)
            scatter!(ax, c.energy_meV, c.intensity; color = :seagreen, markersize = 6,
                     marker = :diamond, label = "corrected (Ei = 4.65)")
        end
        mk = get(model, (B, q), nothing)
        mk === nothing || lines!(ax, mk[1], mk[2]; color = :crimson, linewidth = 2.6,
                                 label = "best-fit model")
        # Raw difference between the two datasets: a first, UNCORRECTED estimate of the
        # Ei = 4.65 magnet background. Resolution is NOT matched -- Ei = 3.32 is sharper --
        # so this over/under-shoots wherever the sample signal is steep. It is the starting
        # point for "healing" the gap, not the answer.
        if hi !== nothing && lo !== nothing
            yl = SV.sv_interp1(lo.energy_meV, lo.intensity, hi.energy_meV)
            d = hi.intensity .- s .* yl
            ok = isfinite.(d) .& (hi.energy_meV .>= 0.35) .& (hi.energy_meV .<= 2.95)
            lines!(ax, hi.energy_meV[ok], d[ok]; color = :purple, linewidth = 1.8,
                   linestyle = :dot, label = "4.65 - 3.32 (resolution NOT matched)")
        end
    end

    vals = Float64[]
    hi === nothing || append!(vals, filter(isfinite, hi.intensity[hi.energy_meV .>= 0.35]))
    lo === nothing || append!(vals, filter(isfinite, s .* lo.intensity[lo.energy_meV .>= 0.35]))
    isempty(vals) || ylims!(ax, -0.08 * maximum(vals), 1.25 * maximum(vals))
    xlims!(ax, 0.0, 4.0)
    row == 1 && col == 1 && axislegend(ax; position = :rt, labelsize = 10)
end

Label(fig[length(fields) + 1, 1:length(cols_q)],
      "Green = the two windows where the background is ANCHORED on min-over-fields data. " *
      "Red = the gap it PCHIP-interpolates across, about 70% of the [0.5, 3.0] meV fit " *
      "window, and where the magnet feature sits. Relative normalisation between the two Ei " *
      "is fitted, so compare positions and shapes only, never absolute heights. Ei = 3.32 " *
      "has the better energy resolution, so a genuinely sharp feature should look sharper " *
      "there.";
      fontsize = 10, color = :grey30)

out = joinpath(FDIR, MODE == "full" ?
    "neutron_1d_two_Ei_with_model_and_background.png" :
    "neutron_1d_two_incident_energies.png")
save(out, fig; px_per_unit = 2)
println("wrote $out")
