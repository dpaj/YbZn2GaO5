#!/usr/bin/env julia
# STAGE 1: build a physics-informed background for the Ei = 3.32 meV data, apply it, and
# overplot the result against the Ei = 4.65 meV corrected cuts. DATA ONLY, ~1 min.
#
#   julia --project=. scripts/background_stage1_ei332.jl
#
# THE ASSUMPTIONS, STATED PLAINLY
#
# Stage 0 showed the field-redundancy test alone only confirms the background above ~2.4 meV
# and below ~0.4 meV at this incident energy -- the middle is unconfirmed because only ONE
# field ever vacates a given energy. This stage therefore adds PHYSICAL priors, and they are
# assumptions rather than measurements, so they are labelled as such throughout:
#
#   A1  There is no known process that produces a peak near 1.5 meV. Anything peaked there is
#       background. (At Ei = 4.65 the corresponding feature sits at ~2.08 meV. A sample
#       excitation cannot move when the incident energy changes, so the shift between the two
#       datasets is itself independent evidence of instrumental origin.)
#   A2  A low-energy upturn is expected and physical in origin -- Bragg tails, incoherent
#       scattering -- so min-over-fields is trustworthy from the low side up to `low_max`.
#   A3  A slow rise toward the kinematic limit is expected -- time-independent pile-up and
#       similar -- so min-over-fields is trustworthy from the high side down to `high_min`.
#   A4  Between `low_max` and `high_min` the background is smooth and is interpolated.
#
# HOW THE COST OF THOSE ASSUMPTIONS IS REPORTED
#
# The window edges and the interpolation form are choices, not measurements, so the script does
# not make one. It builds the background under a GRID of defensible choices and reports the
# envelope. That spread is the "ideological" uncertainty -- the part that comes from picking one
# defensible analysis over another -- and it is propagated into the corrected error bars in
# quadrature with counting statistics, so it travels with the data instead of being an
# after-the-fact caveat.
#
# The output CSV also carries, per point, whether the background there is CONFIRMED by field
# redundancy (Stage 0's statistic) or merely ASSUMED on the grounds above. Those are different
# kinds of evidence and are not merged.

using Printf, Statistics, LinearAlgebra, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const controls = SV.sv_load_controls(REPO_ROOT)
const DIR = SV.sv_repo_path(REPO_ROOT, controls["paths"]["neutron_1d_dir"])
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/background_stage1")
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/background_stage1")
mkpath(FDIR); mkpath(TDIR)

const EI = 3.32
const TBASE = 0.07

# Central choice, and the grid of defensible alternatives whose spread IS the reported
# ideological uncertainty. The by-eye values are low_max = 1.3 and high_min = 1.8.
const LOW_MAX = [1.2, 1.3, 1.4]
const HIGH_MIN = [1.7, 1.8, 1.9]
const INTERP = [:pchip, :linear]
const CENTRAL = (1.3, 1.8, :pchip)

scans = Dict{Tuple{Float64,String},Any}()
for f in sort(readdir(DIR))
    endswith(lowercase(f), ".dat") || continue
    path = joinpath(DIR, f)
    meta = try SV.sv_parse_neutron_1d_filename(path) catch; continue end
    (meta.Ei_meV ≈ EI && meta.temperature_K ≈ TBASE) || continue
    scans[(meta.field_T, meta.qtag)] = SV.sv_load_neutron_raw_scan_1d(path, controls)
end
isempty(scans) && error("No Ei = $EI, T = $TBASE scans in $DIR")

# Ei = 4.65 corrected cuts, for the overplot that validates the whole construction.
cuts465 = Dict{Tuple{Float64,String},Any}()
for c in SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
    cuts465[(c.field_T, c.qtag)] = c
end

qtags = sort(unique(String[k[2] for k in keys(scans)]))
fields = sort(unique(Float64[k[1] for k in keys(scans)]))
qlabel(t) = (p = split(t, "_"); length(p) == 3 ?
             "(" * join(replace.(p, "p" => "."), ", ") * ")" : t)

"Per-energy min and second-min across fields; the gap between them is Stage 0's statistic."
function order_stats(byfield, fs, E)
    lo = fill(NaN, length(E)); lo2 = fill(NaN, length(E))
    for i in eachindex(E)
        vals = Float64[]
        for B in fs
            s = byfield[B]
            y = SV.sv_interp1(s.energy_meV, s.intensity, [E[i]])[1]
            isfinite(y) && push!(vals, y)
        end
        length(vals) >= 2 || continue
        sort!(vals); lo[i] = vals[1]; lo2[i] = vals[2]
    end
    return lo, lo2
end

"Background from min-over-fields anchored outside [low_max, high_min] and interpolated across."
function build_bg(E, lo; low_max, high_min, interp)
    keep = ((E .<= low_max) .| (E .>= high_min)) .& isfinite.(lo)
    count(keep) >= 4 || return fill(NaN, length(E))
    bg, _, _ = SV.sv_make_interpolated_background(E, E[keep], lo[keep];
                   smooth_sigma_meV=0.0, interpolation_kind=interp)
    return Float64.(bg)
end

rows = Vector{Any}()
fig = Figure(size = (640 * length(qtags) + 60, 350 * length(fields) + 130))
Label(fig[0, 1:length(qtags)],
      "Stage 1 -- physics-informed background for Ei = 3.32 meV, applied, and compared with " *
      "the Ei = 4.65 meV corrected cuts";
      fontsize = 16, font = :bold)

for (col, q) in enumerate(qtags)
    byfield = Dict(B => scans[(B, q)] for B in fields if haskey(scans, (B, q)))
    fs = sort(collect(keys(byfield)))
    length(fs) >= 2 || continue
    E = byfield[fs[1]].energy_meV
    lo, lo2 = order_stats(byfield, fs, E)
    confirmed = isfinite.(lo) .& isfinite.(lo2) .& ((lo2 .- lo) ./ max.(abs.(lo), eps()) .< 0.10)

    # Family of backgrounds over the defensible choices; envelope = ideological uncertainty.
    fam = [build_bg(E, lo; low_max=a, high_min=b, interp=c)
           for a in LOW_MAX, b in HIGH_MIN, c in INTERP]
    fam = [f for f in fam if any(isfinite, f)]
    bg_c = build_bg(E, lo; low_max=CENTRAL[1], high_min=CENTRAL[2], interp=CENTRAL[3])
    bg_lo = [minimum(skipmissing([isfinite(f[i]) ? f[i] : missing for f in fam])) for i in eachindex(E)]
    bg_hi = [maximum(skipmissing([isfinite(f[i]) ? f[i] : missing for f in fam])) for i in eachindex(E)]
    ideo = (bg_hi .- bg_lo) ./ 2                      # half-envelope as a 1-sigma-like scale

    @printf("\n%s (%s), %d background variants\n", q, qlabel(q), length(fam))
    for e in (0.8, 1.3, 1.5, 1.8, 2.3)
        i = argmin(abs.(E .- e))
        @printf("  E = %.2f meV: bg = %.3e  ideological +/- %.3e (%.0f%%)   %s\n",
                E[i], bg_c[i], ideo[i], 100 * ideo[i] / max(bg_c[i], eps()),
                confirmed[i] ? "CONFIRMED by field redundancy" : "ASSUMED (A1-A4)")
    end

    for (row, B) in enumerate(fs)
        s = byfield[B]
        ax = Axis(fig[row, col], xlabel = row == length(fs) ? "energy transfer (meV)" : "",
                  ylabel = col == 1 ? "intensity (arb.)" : "",
                  title = "$(qlabel(q))   $(round(Int, B)) T")
        band!(ax, E, bg_lo, bg_hi; color = (:darkorange, 0.22))
        lines!(ax, E, bg_c; color = :darkorange, linewidth = 2.4, linestyle = :dashdot,
               label = "background (band = assumption spread)")
        scatter!(ax, s.energy_meV, s.intensity; color = (:grey40, 0.75), markersize = 5,
                 label = "raw Ei = 3.32")
        # Corrected, with the ideological spread propagated into the error bar.
        corr = s.intensity .- bg_c
        etot = sqrt.(s.error .^ 2 .+ ideo .^ 2)
        errorbars!(ax, E, corr, etot; color = (:black, 0.4), whiskerwidth = 3)
        scatter!(ax, E, corr; color = :black, markersize = 6, label = "corrected Ei = 3.32")
        # The validation: does it agree with the independently corrected 4.65 data?
        c465 = get(cuts465, (B, q), nothing)
        if c465 !== nothing
            y = SV.sv_interp1(c465.energy_meV, c465.intensity, E)
            ok = isfinite.(y) .& isfinite.(corr) .& (E .>= 0.5) .& (E .<= 2.4)
            sc = count(ok) > 4 ? max(0.0, dot(y[ok], corr[ok]) / max(eps(), dot(y[ok], y[ok]))) : 1.0
            lines!(ax, c465.energy_meV, sc .* c465.intensity; color = :dodgerblue,
                   linewidth = 2.4, label = @sprintf("corrected Ei = 4.65 (x%.3g)", sc))
        end
        vspan!(ax, CENTRAL[1], CENTRAL[2]; color = (:indianred, 0.08))
        xlims!(ax, 0.0, 3.0)
        vv = filter(isfinite, vcat(corr, s.intensity[s.energy_meV .>= 0.5]))
        isempty(vv) || ylims!(ax, min(-0.0004, 1.3 * minimum(vv)), 1.2 * maximum(vv))
        row == 1 && col == 1 && axislegend(ax; position = :rt, labelsize = 9)

        for i in eachindex(E)
            push!(rows, (q, B, E[i], s.intensity[i], bg_c[i], bg_lo[i], bg_hi[i],
                         corr[i], s.error[i], etot[i], confirmed[i]))
        end
    end
end

open(joinpath(TDIR, "ei3p32_background_applied.csv"), "w") do io
    println(io, "qtag,field_T,energy_meV,raw,bg_central,bg_envelope_lo,bg_envelope_hi," *
                "corrected,err_counting,err_total,background_confirmed_by_redundancy")
    for r in rows
        @printf(io, "%s,%.1f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n", r...)
    end
end

Label(fig[length(fields) + 1, 1:length(qtags)],
      "Red shading is the interpolated region, where the background is ASSUMED rather than " *
      "confirmed. Orange band is the spread over 18 defensible choices of window edges and " *
      "interpolation form -- the ideological uncertainty -- and it is propagated into the " *
      "corrected error bars in quadrature with counting statistics. The blue curve is the " *
      "independently corrected Ei = 4.65 data: the two constructions share no assumptions, so " *
      "agreement between them is a real test.";
      fontsize = 10, color = :grey30)

out = joinpath(FDIR, "background_stage1_ei3p32.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
println("wrote " * joinpath(TDIR, "ei3p32_background_applied.csv"))
