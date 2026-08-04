#!/usr/bin/env julia
# STAGE 1 for the PRODUCTION Ei = 4.65 meV data: build the background under a family of defensible
# choices, and report the envelope as a per-point uncertainty. DATA ONLY, ~1 min.
#
#   julia --project=. scripts/background_stage1_ei465.jl
#
# This is the input the neutron objective needs in order to weight energies by how well the
# background is known there, instead of by counting statistics alone. Measured motivation: on the
# two 9 T dispersive cuts, 35-41% of chi2 comes from the 1.8-2.4 meV band where their modes are NOT
# -- i.e. from an instrumental artefact -- and that is what pulled gzz toward 3.70.
#
# THE DESIGN PROBLEM, AND WHY VARYING WINDOW EDGES ALONE WOULD LIE
#
# The production construction anchors min-over-fields on [0, 0.75] meV and above 2.5 meV and
# PCHIP-interpolates the gap between. The ~2.08 meV magnet feature sits inside that gap. So if the
# variant family only moves the window edges and the interpolation form, EVERY member smooths
# straight across the feature, they all agree closely, and the envelope comes out small -- while
# being wrong in the same direction. A band that says "we agree" when what we agree on is wrong is
# worse than no band.
#
# So the family deliberately BRACKETS the feature instead:
#
#   FAMILY A (smooth)  -- interpolate across the gap, as production does. The magnet feature is
#                         left IN the data, so this UNDER-subtracts: a lower bound on the
#                         background.
#   FAMILY B (min-fed) -- use min-over-fields inside the gap too. The magnet background is
#                         instrumental and therefore present at EVERY field, so min-over-fields
#                         does not remove it -- it is only excluded because the production windows
#                         skip the gap. But in the gap the minimum still contains signal from
#                         whichever field has least there, so this OVER-subtracts: an upper bound.
#
# The truth lies between, and the envelope over both families is an honest statement of that. This
# is the same principle as the Ei = 3.32 stage, with one addition: there the assumption cost came
# only from window placement, here it is dominated by whether the feature is modelled at all.
#
# `(0,1,0)` is treated slightly differently and the asymmetry is deliberate. Its 0 T scan shows no
# 2.08 meV feature, so for that cut the two families should nearly coincide -- and if they do, that
# is independent confirmation that the cut is clean, which is what makes it the trustworthy gzz
# determination. If they do NOT coincide, that conclusion needs revisiting.

using Printf, Statistics, LinearAlgebra, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const controls = SV.sv_load_controls(REPO_ROOT)
const DIR = SV.sv_repo_path(REPO_ROOT, controls["paths"]["neutron_1d_dir"])
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/background_stage1_ei465")
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/background_stage1_ei465")
mkpath(FDIR); mkpath(TDIR)

const EI = 4.65
const TBASE = 0.07
# Centred on the production windows: [0, 0.75] and above 2.5 meV.
const LOW_MAX = [0.65, 0.75, 0.85]
const HIGH_MIN = [2.3, 2.5, 2.7]
const INTERP = [:pchip, :linear]
const CENTRAL = (0.75, 2.5, :pchip)

scans = Dict{Tuple{Float64,String},Any}()
for f in sort(readdir(DIR))
    endswith(lowercase(f), ".dat") || continue
    path = joinpath(DIR, f)
    meta = try SV.sv_parse_neutron_1d_filename(path) catch; continue end
    (meta.Ei_meV ≈ EI && meta.temperature_K ≈ TBASE) || continue
    scans[(meta.field_T, meta.qtag)] = SV.sv_load_neutron_raw_scan_1d(path, controls)
end
isempty(scans) && error("No Ei = $EI, T = $TBASE scans in $DIR")
qtags = sort(unique(String[k[2] for k in keys(scans)]))
fields = sort(unique(Float64[k[1] for k in keys(scans)]))
@printf("Ei = %.2f meV: %d qtags, %d fields %s\n", EI, length(qtags), length(fields), string(fields))

qlabel(t) = (p = split(t, "_"); length(p) == 3 ?
             "(" * join(replace.(p, "p" => "."), ", ") * ")" : t)

"Per-energy min and second-min across fields. The gap between them is the Stage 0 statistic."
function order_stats(byfield, fs, E)
    lo = fill(NaN, length(E)); lo2 = fill(NaN, length(E))
    for i in eachindex(E)
        v = Float64[]
        for B in fs
            y = SV.sv_interp1(byfield[B].energy_meV, byfield[B].intensity, [E[i]])[1]
            isfinite(y) && push!(v, y)
        end
        length(v) >= 2 || continue
        sort!(v); lo[i] = v[1]; lo2[i] = v[2]
    end
    return lo, lo2
end

"""
Background under one choice. `family = :smooth` anchors only outside [low_max, high_min] and
interpolates the gap (production behaviour, a LOWER bound). `family = :minfed` also feeds the
min-over-fields values inside the gap to the interpolator (an UPPER bound).
"""
function build_bg(E, lo; low_max, high_min, interp, family)
    keep = isfinite.(lo) .& (family === :minfed ? trues(length(E)) :
                             ((E .<= low_max) .| (E .>= high_min)))
    count(keep) >= 4 || return fill(NaN, length(E))
    bg, _, _ = SV.sv_make_interpolated_background(E, E[keep], lo[keep];
                   smooth_sigma_meV=0.0, interpolation_kind=interp)
    return Float64.(bg)
end

rows = Vector{Any}()
fig = Figure(size = (600 * length(qtags) + 60, 340 * length(fields) + 150))
Label(fig[0, 1:length(qtags)],
      "Stage 1, Ei = 4.65 meV: background under 36 defensible choices. The band BRACKETS whether " *
      "the 2.08 meV magnet feature is modelled.";
      fontsize = 15, font = :bold)

for (col, q) in enumerate(qtags)
    byfield = Dict(B => scans[(B, q)] for B in fields if haskey(scans, (B, q)))
    fs = sort(collect(keys(byfield)))
    length(fs) >= 2 || continue
    E = byfield[fs[1]].energy_meV
    lo, lo2 = order_stats(byfield, fs, E)
    confirmed = isfinite.(lo) .& isfinite.(lo2) .&
                ((lo2 .- lo) ./ max.(abs.(lo), eps()) .< 0.10)

    famA = [build_bg(E, lo; low_max=a, high_min=b, interp=c, family=:smooth)
            for a in LOW_MAX, b in HIGH_MIN, c in INTERP]
    famB = [build_bg(E, lo; low_max=a, high_min=b, interp=c, family=:minfed)
            for a in LOW_MAX, b in HIGH_MIN, c in INTERP]
    all_v = [f for f in vcat(vec(famA), vec(famB)) if any(isfinite, f)]
    bg_c = build_bg(E, lo; low_max=CENTRAL[1], high_min=CENTRAL[2], interp=CENTRAL[3],
                    family=:smooth)
    finite_at(i) = [f[i] for f in all_v if isfinite(f[i])]
    bg_lo = [isempty(finite_at(i)) ? NaN : minimum(finite_at(i)) for i in eachindex(E)]
    bg_hi = [isempty(finite_at(i)) ? NaN : maximum(finite_at(i)) for i in eachindex(E)]
    ideo = (bg_hi .- bg_lo) ./ 2

    # How much of the envelope comes from the FAMILY choice rather than the window placement?
    aA = [f for f in vec(famA) if any(isfinite, f)]
    aB = [f for f in vec(famB) if any(isfinite, f)]
    spreadA = [isempty([f[i] for f in aA if isfinite(f[i])]) ? NaN :
               maximum([f[i] for f in aA if isfinite(f[i])]) -
               minimum([f[i] for f in aA if isfinite(f[i])]) for i in eachindex(E)]
    famgap = [isempty([f[i] for f in aB if isfinite(f[i])]) || isnan(bg_c[i]) ? NaN :
              median([f[i] for f in aB if isfinite(f[i])]) - bg_c[i] for i in eachindex(E)]

    @printf("\n%s (%s), %d variants\n", q, qlabel(q), length(all_v))
    for e in (0.9, 1.5, 2.08, 2.4)
        i = argmin(abs.(E .- e))
        @printf("  E = %.2f: bg = %.3e  envelope +/- %.3e (%.0f%%)  window-only %.3e  family gap %.3e  %s\n",
                E[i], bg_c[i], ideo[i], 100*ideo[i]/max(bg_c[i], eps()),
                spreadA[i]/2, famgap[i],
                confirmed[i] ? "CONFIRMED" : "assumed")
    end
    fw = filter(isfinite, spreadA[(E .>= 0.75) .& (E .<= 2.5)]) ./ 2
    fg = filter(isfinite, abs.(famgap[(E .>= 0.75) .& (E .<= 2.5)]))
    isempty(fw) || @printf("  across the gap: window-placement spread %.3e, FAMILY gap %.3e (%.1fx larger)\n",
                           median(fw), median(fg), median(fg)/max(median(fw), eps()))

    for (row, B) in enumerate(fs)
        s = byfield[B]
        ax = Axis(fig[row, col], xlabel = row == length(fs) ? "energy transfer (meV)" : "",
                  ylabel = col == 1 ? "intensity (arb.)" : "",
                  title = "$(qlabel(q))   $(round(Int, B)) T")
        vspan!(ax, CENTRAL[1], CENTRAL[2]; color = (:indianred, 0.07))
        band!(ax, E, bg_lo, bg_hi; color = (:darkorange, 0.25))
        lines!(ax, E, bg_c; color = :darkorange, linewidth = 2.4, linestyle = :dashdot,
               label = "background, central choice")
        scatter!(ax, s.energy_meV, s.intensity; color = (:grey40, 0.7), markersize = 5,
                 label = "raw")
        corr = s.intensity .- bg_c
        etot = sqrt.(s.error .^ 2 .+ ideo .^ 2)
        errorbars!(ax, E, corr, etot; color = (:black, 0.35), whiskerwidth = 3)
        scatter!(ax, E, corr; color = :black, markersize = 5, label = "corrected")
        xlims!(ax, 0.0, 3.4)
        vv = filter(isfinite, vcat(corr, s.intensity[s.energy_meV .>= 0.5]))
        isempty(vv) || ylims!(ax, min(-0.0004, 1.3*minimum(vv)), 1.15*maximum(vv))
        row == 1 && col == 1 && axislegend(ax; position = :rt, labelsize = 9)

        for i in eachindex(E)
            push!(rows, (q, B, E[i], s.intensity[i], bg_c[i], bg_lo[i], bg_hi[i],
                         ideo[i], corr[i], s.error[i], etot[i], confirmed[i]))
        end
    end
end

open(joinpath(TDIR, "ei4p65_background_envelope.csv"), "w") do io
    println(io, "qtag,field_T,energy_meV,raw,bg_central,bg_envelope_lo,bg_envelope_hi," *
                "bg_sigma,corrected,err_counting,err_total,background_confirmed_by_redundancy")
    for r in rows
        @printf(io, "%s,%.1f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n", r...)
    end
end

Label(fig[length(fields) + 1, 1:length(qtags)],
      "Orange band = envelope over 36 choices: 3 low edges x 3 high edges x 2 interpolation forms x " *
      "2 FAMILIES. The families bracket the 2.08 meV magnet feature -- smooth interpolation leaves " *
      "it in the data (under-subtracts), min-over-fields inside the gap removes it along with some " *
      "signal (over-subtracts). Varying window edges alone would have produced a small band in " *
      "which every member was wrong the same way. `bg_sigma` is half the envelope and is what the " *
      "objective should add to counting statistics in quadrature.";
      fontsize = 10, color = :grey30)

out = joinpath(FDIR, "background_stage1_ei465.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
println("wrote " * joinpath(TDIR, "ei4p65_background_envelope.csv"))
