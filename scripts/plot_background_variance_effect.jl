#!/usr/bin/env julia
# What the background-variance term does to the gzz determination. READS CSV ONLY -- no Sunny, no
# recompute, ~20 s. Generate the input with scripts/check_background_variance_effect.jl (~1 h).
#
#   julia --project=. scripts/plot_background_variance_effect.jl
#
# THE ACCEPTANCE TEST FAILED ITS OWN PREDICTION, AND THAT IS THE RESULT
#
# The prediction was specific and falsifiable: the two 9 T DISPERSIVE cuts should move OFF gzz =
# 3.70 once the 1.8-2.4 meV band gets large error bars, because 35-41% of their chi2 came from that
# band and the 2.08 meV magnet background lives there. THEY DO NOT MOVE AT ALL. Both still prefer
# 3.70. Only the two 14 T cuts shift, and only by -0.10.
#
# It is not that the term is inert or mis-targeted. It is strongly selective and it does a great
# deal to chi2:
#   * the envelope sigma for (0.33,0.33,0) is 5.8e-4 at 2.02 meV against 5.0e-5 at 1.03 meV where
#     its mode actually sits -- a factor 11.5 in sigma, 133 in VARIANCE, aimed exactly at the
#     artefact and away from the signal;
#   * in-window variance inflation is 5.2x and the six-cut chi2_red falls 19.74 -> 5.58.
#
# WHY THE MINIMUM DOES NOT MOVE, and this is the part worth keeping. Down-weighting the band scales
# each cut's whole chi2(gzz) CURVE down -- level to 17-35% and curvature to 7-48% -- rather than
# subtracting a gzz-independent constant from it. A constant would shift nothing; a uniform scaling
# also shifts nothing. Either way the argmin survives.
#
# SO AN EARLIER INFERENCE HERE WAS WRONG. "35-41% of chi2 comes from the artefact band, and that is
# the mechanism that pulled gzz toward 3.70" does not follow, and is now refuted directly: removing
# that band's influence by 133x in variance leaves the preferred gzz exactly where it was. The chi2
# FRACTION sitting in a band is the wrong diagnostic for what DRIVES a parameter. What matters is the
# band's contribution to the DERIVATIVE of chi2 with respect to the parameter, and by that measure
# the artefact band was never the driver.
#
# CONSEQUENCE, and it cuts in our favour: the dispersive cuts' preference for gzz = 3.70 is NOT a
# background artefact. The tension with (0,1,0) is therefore real and needs a different explanation
# -- exchange anisotropy, a genuinely q-dependent systematic, or model error. It is no longer
# attributable to the background.
#
# TWO CAVEATS THAT LIMIT THIS RUN. (0,1,0) at 9 T prefers gzz = 3.30, which is the LOW EDGE of the
# scanned range, so its preferred value is not determined -- only bounded above. A wider scan is
# needed before quoting it. And with the variance on, (0,1,0) at 14 T is nearly flat (span 0.31 over
# the whole range), i.e. it no longer constrains gzz at all; that is arguably CORRECT, since its
# mode exits the 3.0 meV fit window in this range and it should not be voting.

using Printf, Statistics, DelimitedFiles, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const TDIR = SV.sv_repo_path(REPO_ROOT,
    "results/feature_tables/sunny_validation/background_variance_effect")
const FDIR = SV.sv_repo_path(REPO_ROOT,
    "results/figures/sunny_validation/background_variance_effect")
mkpath(FDIR)
const CSV = joinpath(TDIR, "gzz_scan_variance_on_off.csv")
isfile(CSV) || error("No scan table. Run scripts/check_background_variance_effect.jl first: $CSV")

raw, hdr = readdlm(CSV, ','; header=true)
ix = Dict(strip(String(h)) => i for (i, h) in enumerate(vec(hdr)))
col(n) = raw[:, ix[n]]
wt   = String.(strip.(String.(col("weighting"))))
gzz  = Float64.(col("gzz"))
tot  = Float64.(col("chi2_red_total"))
infl = Float64.(col("inflation"))
B    = Float64.(col("field_T"))
qt   = String.(strip.(String.(col("qtag"))))
cut  = Float64.(col("chi2_red_cut"))

G = sort(unique(gzz))
cuts = unique(collect(zip(qt, B)))
sixcut(mode) = [only(unique(tot[(wt .== mode) .& (gzz .≈ g)])) for g in G]
percut(mode, q, b) = [only(cut[(wt .== mode) .& (gzz .≈ g) .& (qt .== q) .& (B .≈ b)]) for g in G]

fig = Figure(size = (1560, 1000))
Label(fig[0, 1:3],
      "Does propagating the background uncertainty change which gzz the data prefer?  " *
      "NO — and that refutes the mechanism it was built on.";
      fontsize = 16, font = :bold)

for (k, (q, b)) in enumerate(sort(cuts, by = x -> (x[2], x[1])))
    r, c = divrem(k - 1, 3)
    off, on = percut("off", q, b), percut("on", q, b)
    ax = Axis(fig[1 + r, 1 + c], xlabel = "gzz",
              ylabel = c == 0 ? "per-cut chi2_red − own minimum" : "",
              title = @sprintf("%s  @ %.0f T", q, b))
    lines!(ax, G, off .- minimum(off); color = :grey35, linewidth = 2.8, label = "variance OFF")
    scatter!(ax, G, off .- minimum(off); color = :grey35, markersize = 9)
    lines!(ax, G, on .- minimum(on); color = :crimson, linewidth = 2.8, label = "variance ON")
    scatter!(ax, G, on .- minimum(on); color = :crimson, markersize = 9)
    goff, gon = G[argmin(off)], G[argmin(on)]
    vlines!(ax, [goff]; color = (:grey35, 0.55), linestyle = :dash)
    vlines!(ax, [gon]; color = (:crimson, 0.55), linestyle = :dot)
    edge = (gon == first(G) || gon == last(G)) ? "  EDGE" : ""
    text!(ax, 0.03, 0.94; space = :relative, fontsize = 10,
          color = goff == gon ? :grey20 : :darkorange4,
          text = goff == gon ? @sprintf("both prefer %.2f%s", gon, edge) :
                              @sprintf("%.2f → %.2f  MOVED%s", goff, gon, edge))
    text!(ax, 0.03, 0.80; space = :relative, fontsize = 9, color = :grey35,
          text = @sprintf("level %.0f%% of OFF, curvature %.0f%%",
                          100*minimum(on)/minimum(off),
                          100*(maximum(on)-minimum(on))/(maximum(off)-minimum(off))))
    k == 1 && axislegend(ax; position = :rt, labelsize = 9)
end

ax = Axis(fig[3, 1:2], xlabel = "gzz", ylabel = "six-cut chi2_red",
          title = "Six-cut objective — the level drops 3.5x, the minimum stays at 3.60",
          yscale = log10)
so, sn = sixcut("off"), sixcut("on")
lines!(ax, G, so; color = :grey35, linewidth = 3, label = "variance OFF")
scatter!(ax, G, so; color = :grey35, markersize = 10)
lines!(ax, G, sn; color = :crimson, linewidth = 3, label = "variance ON")
scatter!(ax, G, sn; color = :crimson, markersize = 10)
vlines!(ax, [G[argmin(so)]]; color = (:grey35, 0.5), linestyle = :dash)
vlines!(ax, [G[argmin(sn)]]; color = (:crimson, 0.5), linestyle = :dot)
axislegend(ax; position = :rt, labelsize = 10)

# NOTE: the @sprintf format string stays a SINGLE LITERAL and the prose is concatenated OUTSIDE the
# macro call -- see the @printf gotcha in CLAUDE.md.
Label(fig[3, 3],
      @sprintf("Median in-window variance inflation: %.2fx", only(unique(infl[wt .== "on"]))) *
      "\n\nThe term is strongly SELECTIVE —\nfor (0.33,0.33,0) the envelope sigma\n" *
      "is 5.8e-4 at 2.02 meV against\n5.0e-5 at 1.03 meV where its mode\n" *
      "is: 11.5x in sigma, 133x in variance.\n\n" *
      "Yet it scales each chi2(gzz) curve\nrather than subtracting a\n" *
      "gzz-dependent piece, so the argmin\nsurvives. The chi2 FRACTION in a\n" *
      "band does not tell you what drives\na parameter — its contribution to\n" *
      "the chi2 DERIVATIVE does.",
      fontsize = 10, color = :grey25, halign = :left, justification = :left, tellwidth = false)

out = joinpath(FDIR, "background_variance_effect.png")
save(out, fig; px_per_unit = 2)

println("six-cut minimum: OFF gzz=", G[argmin(so)], "  ON gzz=", G[argmin(sn)])
@printf("per-cut: %d of %d cuts changed preference\n",
        count(((q, b),) -> G[argmin(percut("off", q, b))] != G[argmin(percut("on", q, b))], cuts),
        length(cuts))
for (q, b) in sort(cuts, by = x -> (x[2], x[1]))
    off, on = percut("off", q, b), percut("on", q, b)
    @printf("  %-14s %5.0f T   OFF %.2f -> ON %.2f   level %3.0f%%  curvature %3.0f%%%s\n",
            q, b, G[argmin(off)], G[argmin(on)],
            100*minimum(on)/minimum(off),
            100*(maximum(on)-minimum(on))/(maximum(off)-minimum(off)),
            G[argmin(on)] in (first(G), last(G)) ? "   <- AT SCAN EDGE, not determined" : "")
end
println("\nwrote $out")
