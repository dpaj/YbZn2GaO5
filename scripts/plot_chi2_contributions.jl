#!/usr/bin/env julia
# Which energies actually drive the fit? Per-point chi2 contribution, (model - data)^2 / sigma^2,
# for each cut and each parameter set. Reads the saved model curves -- NO Sunny, NO KPM, ~1 min.
#
#   julia --project=. scripts/plot_chi2_contributions.jl
#
# WHY THIS IS THE PLOT TO LOOK AT BEFORE CHANGING THE OBJECTIVE
#
# chi2_red is a single number and it hides where it came from. The fit does not weight energies
# by how much we trust them -- it weights by 1/sigma^2 with sigma from counting statistics alone,
# so a well-counted point in a badly-known background region carries FULL weight and a
# poorly-counted point in a well-known region carries little. That is backwards relative to what
# we actually know.
#
# The specific worry: the ~2.08 meV magnet background sits inside the PCHIP-interpolated gap that
# the min-over-fields construction cannot reach, so it survives into the corrected data. If that
# region contributes a large share of the total chi2, then the fit is substantially being driven
# by an instrumental artefact -- which is the mechanism that pulled gzz toward 3.70. This
# quantifies that share instead of asserting it.
#
# The cumulative curve is the useful half: it answers "what fraction of this cut's chi2 comes from
# below energy E", so a step in it locates the energies doing the work.

using Printf, Statistics, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const CSV = SV.sv_repo_path(REPO_ROOT,
    "results/feature_tables/sunny_validation/neutron_parameter_sets/model_curves.csv")
isfile(CSV) || error("No $CSV; run scripts/plot_neutron_parameter_sets.jl first.")
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/chi2_contributions")
mkpath(FDIR)

# The energy band containing the magnet background, from the two-incident-energy comparison:
# the feature sits at ~2.08 meV at Ei = 4.65 with a width of order the resolution.
const BG_BAND = (1.80, 2.40)
# Approximate magnon position per (qtag, field) at the fitted parameters, from the overplots.
# Used only for annotation, so the band's chi2 share can be interpreted rather than just read.
const MODE_POS = Dict(
    ("0_1_0", 9.0) => 1.90, ("0_1_0", 14.0) => 3.05,
    ("0p33_0p33_0", 9.0) => 1.05, ("0p33_0p33_0", 14.0) => 2.05,
    ("0p5_0_0", 9.0) => 1.10, ("0p5_0_0", 14.0) => 2.00,
)
const WINDOW = (0.5, 3.0)

lines_ = filter(!isempty, strip.(readlines(CSV)))
hdr = String.(split(lines_[1], ','))
ix = Dict(h => i for (i, h) in enumerate(hdr))
rows = [String.(split(l, ',')) for l in lines_[2:end]]
rows = [r for r in rows if length(r) == length(hdr)]
num(r, k) = something(tryparse(Float64, r[ix[k]]), NaN)

sets = unique(String[r[ix["set"]] for r in rows])
qtags = unique(String[r[ix["qtag"]] for r in rows])
fields = sort(unique(Float64[num(r, "field_T") for r in rows]))
want = let e = get(ENV, "YZGO_CHI2_SETS", "")
    isempty(e) ? sets : [s for s in sets if s in String.(strip.(split(e, ',')))]
end

qlabel(t) = (p = split(t, "_"); length(p) == 3 ?
             "(" * join(replace.(p, "p" => "."), ", ") * ")" : t)

"Per-point chi2 contribution inside the fit window, plus the cumulative fraction."
function contributions(set, q, B)
    sel = [r for r in rows if r[ix["set"]] == set && r[ix["qtag"]] == q &&
                              num(r, "field_T") ≈ B]
    E = num.(sel, "energy_meV"); y = num.(sel, "I_exp")
    e = num.(sel, "Ierr_exp");   m = num.(sel, "I_model_scaled")
    p = sortperm(E); E, y, e, m = E[p], y[p], e[p], m[p]
    inwin = (E .>= WINDOW[1]) .& (E .<= WINDOW[2]) .& isfinite.(y) .& isfinite.(m)
    c = fill(NaN, length(E))
    for i in eachindex(E)
        inwin[i] || continue
        w = (isfinite(e[i]) && e[i] > 0) ? 1 / e[i]^2 : 1.0
        c[i] = w * (m[i] - y[i])^2
    end
    tot = sum(filter(isfinite, c))
    cum = cumsum([isfinite(v) ? v : 0.0 for v in c]) ./ max(tot, eps())
    return (; E, c, cum, tot, inwin)
end

@printf("%-12s %-14s %6s %10s %12s %10s\n", "set", "qtag", "field", "chi2_sum",
        "in 1.8-2.4", "share")
shares = Dict{String,Vector{Float64}}()
for set in want, q in qtags, B in fields
    r = contributions(set, q, B)
    isfinite(r.tot) && r.tot > 0 || continue
    band = sum(v for (E, v) in zip(r.E, r.c)
               if isfinite(v) && BG_BAND[1] <= E <= BG_BAND[2]; init=0.0)
    @printf("%-12s %-14s %5.0f T %10.1f %12.1f %9.0f%%\n", set, qlabel(q), B, r.tot,
            band, 100 * band / r.tot)
    push!(get!(shares, set, Float64[]), 100 * band / r.tot)
end
println()
for set in want
    haskey(shares, set) || continue
    v = shares[set]
    @printf("%-12s  the %.1f-%.1f meV band carries %.0f%% of chi2 on average (range %.0f-%.0f%%), from %.0f%% of the window width
",
            set, BG_BAND[1], BG_BAND[2], mean(v), minimum(v), maximum(v),
            100 * (BG_BAND[2] - BG_BAND[1]) / (WINDOW[2] - WINDOW[1]))
end

cols = [:crimson, :dodgerblue, :seagreen, :darkorange]
fig = Figure(size = (540 * length(qtags) + 60, 340 * length(fields) + 120))
Label(fig[0, 1:length(qtags)],
      "Where the fit's chi2 actually comes from. Bars = per-point (model - data)^2 / sigma^2; " *
      "line = cumulative fraction. Shaded = the magnet-background band.";
      fontsize = 15, font = :bold)

for (row, B) in enumerate(fields), (col, q) in enumerate(qtags)
    ax = Axis(fig[row, col], xlabel = row == length(fields) ? "energy transfer (meV)" : "",
              ylabel = col == 1 ? "chi2 contribution" : "",
              title = "$(qlabel(q))   $(round(Int, B)) T", yscale = log10)
    ax2 = Axis(fig[row, col], yaxisposition = :right, ylabel = col == length(qtags) ?
               "cumulative fraction" : "", ygridvisible = false)
    hidespines!(ax2); hidexdecorations!(ax2)
    linkxaxes!(ax, ax2)
    vspan!(ax, BG_BAND[1], BG_BAND[2]; color = (:indianred, 0.13))
    # Where the real mode sits, so the band's chi2 share can be read correctly: a large share is
    # EXPECTED and desirable where the magnon is inside the band, and a warning sign only where
    # the magnon is elsewhere and the band therefore holds artefact rather than signal.
    mode_E = get(MODE_POS, (q, B), NaN)
    isfinite(mode_E) && vlines!(ax, [mode_E]; color = :black, linewidth = 2.5,
                               linestyle = :dot)
    for (k, set) in enumerate(want)
        r = contributions(set, q, B)
        isfinite(r.tot) && r.tot > 0 || continue
        cc = [isfinite(v) && v > 0 ? v : NaN for v in r.c]
        scatter!(ax, r.E, cc; color = (cols[mod1(k, length(cols))], 0.85), markersize = 6,
                 label = @sprintf("%s (chi2 = %.0f)", set, r.tot))
        lines!(ax2, r.E, r.cum; color = cols[mod1(k, length(cols))], linewidth = 2.4,
               linestyle = :dash)
    end
    allc = Float64[]
    for set in want
        r = contributions(set, q, B)
        append!(allc, filter(v -> isfinite(v) && v > 0, r.c))
    end
    if !isempty(allc)
        hi = maximum(allc)
        ylims!(ax, max(1e-8, hi * 1e-5), hi * 3)
    end
    xlims!(ax, 0.3, 3.2); xlims!(ax2, 0.3, 3.2); ylims!(ax2, 0, 1.05)
    row == 1 && col == 1 && axislegend(ax; position = :lt, labelsize = 9)
end

Label(fig[length(fields) + 1, 1:length(qtags)],
      "Dotted vertical line = where the magnon sits. READ THE SHARE AGAINST IT: a large share " *
      "inside the shaded band is expected and desirable where the mode is IN the band, and a " *
      "warning only where the mode is elsewhere so the band holds artefact rather than signal. " *
      "The fit weights by 1/sigma^2 from counting statistics alone, so a well-counted point in a " *
      "badly-known background region carries full weight regardless.";
      fontsize = 10, color = :grey30)

out = joinpath(FDIR, "chi2_contributions.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
