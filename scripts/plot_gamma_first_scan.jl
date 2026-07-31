#!/usr/bin/env julia
# Diagnostic overplots for the Gamma-first parameter scan. READS CSV ONLY -- no KPM, no
# Sunny calls, so it is free to run while a scan is using the box.
#
#   julia --project=. scripts/plot_gamma_first_scan.jl
#
# Panels:
#   1  Gamma chi2 surface in (gzz, sigma_gzz), two runs merged to span gzz 2.85-4.05
#   2  Gamma marginals -- are gzz and sigma_gzz interior minima, or grid artefacts?
#   3  K/M chi2 surface in (J1, sigma_J)
#   4  K/M chi2 vs sigma_J, one line per J1. THE key panel: sigma_J lowers chi2 at EVERY
#      J1, i.e. it absorbs misfit rather than trading off against the bandwidth.
#   5  Stage-4 joint refinement trajectory on all six cuts
#   6  Per-cut chi2 -- the zone centre fits ~25-50x better than the dispersive cuts

using Printf, Statistics, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const TDIR = joinpath(REPO_ROOT, "results", "feature_tables", "sunny_validation",
                      "gamma_first_scan")
const FDIR = joinpath(REPO_ROOT, "results", "figures", "sunny_validation",
                      "gamma_first_scan")
mkpath(FDIR)

function readcsv(path)
    isfile(path) || return (String[], Dict{String,String}[])
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && return (String[], Dict{String,String}[])
    hdr = String.(split(lines[1], ','))
    rows = Dict{String,String}[]
    for l in lines[2:end]
        f = String.(split(l, ','))
        length(f) == length(hdr) || continue
        push!(rows, Dict(zip(hdr, f)))
    end
    return (hdr, rows)
end
num(r, k) = (v = get(r, k, ""); isempty(v) ? NaN : something(tryparse(Float64, v), NaN))
finite(v) = filter(isfinite, v)

"Grid a set of rows onto a matrix of chi2, NaN where a point was never evaluated."
function surface(rows, xkey, ykey; vkey="chi2_red")
    xs = sort(unique(finite([num(r, xkey) for r in rows])))
    ys = sort(unique(finite([num(r, ykey) for r in rows])))
    Z = fill(NaN, length(xs), length(ys))
    for r in rows
        x, y, v = num(r, xkey), num(r, ykey), num(r, vkey)
        (isfinite(x) && isfinite(y) && isfinite(v)) || continue
        i = findfirst(≈(x), xs); j = findfirst(≈(y), ys)
        (i === nothing || j === nothing) && continue
        # Keep the better value if two runs overlap on the same point.
        Z[i, j] = isnan(Z[i, j]) ? v : min(Z[i, j], v)
    end
    return (xs, ys, Z)
end

argmin_grid(xs, ys, Z) = begin
    best = (Inf, NaN, NaN)
    for i in eachindex(xs), j in eachindex(ys)
        isfinite(Z[i, j]) && Z[i, j] < best[1] && (best = (Z[i, j], xs[i], ys[j]))
    end
    best
end

_, s1a = readcsv(joinpath(TDIR, "run1_suspended", "stage1_gamma_gzz_sigma_gzz.csv"))
_, s1b = readcsv(joinpath(TDIR, "stage1_gamma_gzz_sigma_gzz.csv"))
s2hdr, s2 = readcsv(joinpath(TDIR, "stage2_disp_J1_sigma_J.csv"))
_, s3 = readcsv(joinpath(TDIR, "stage3_gamma_recheck.csv"))
_, s4 = readcsv(joinpath(TDIR, "stage4_joint_refinement.csv"))

fig = Figure(size = (1680, 980))
Label(fig[0, 1:3],
      "YbZn2GaO5 -- Gamma-first parameter scan on the 1D neutron cuts, 36x36x1, 81 q, " *
      "4 realizations (common random numbers)";
      fontsize = 17, font = :bold)

# ---------------------------------------------------------------- 1  Gamma surface
gam = vcat(s1a, s1b)
gx, gy, GZ = surface(gam, "gzz", "sigma_gzz")
ax1 = Axis(fig[1, 1], xlabel = "gzz", ylabel = "sigma_gzz",
           title = "(0,1,0) only: J cancels at q=0, so this is a pure g probe")
hm1 = heatmap!(ax1, gx, gy, log10.(GZ); colormap = :viridis)
Colorbar(fig[1, 1][1, 2], hm1, label = "log10 chi2_red")
gbest = argmin_grid(gx, gy, GZ)
scatter!(ax1, [gbest[2]], [gbest[3]]; marker = :star5, markersize = 26,
         color = :white, strokecolor = :black, strokewidth = 1.5)
scatter!(ax1, [3.8], [0.8]; marker = :xcross, markersize = 18, color = :red,
         strokecolor = :white, strokewidth = 1)
vlines!(ax1, [3.44]; color = (:orange, 0.9), linestyle = :dash, linewidth = 2)
text!(ax1, 3.46, 0.05; text = "published\ngzz=3.44", color = :orange, fontsize = 11)
text!(ax1, 3.72, 0.86; text = "by-eye", color = :red, fontsize = 11)
@printf("Gamma minimum: gzz = %.2f, sigma_gzz = %.2f, chi2_red = %.2f\n",
        gbest[2], gbest[3], gbest[1])

# ---------------------------------------------------------------- 2  Gamma marginals
ax2 = Axis(fig[1, 2], xlabel = "gzz", ylabel = "chi2_red on (0,1,0)",
           title = "gzz has an INTERIOR minimum at every sigma_gzz", yscale = log10)
for (j, sg) in enumerate(gy)
    v = [GZ[i, j] for i in eachindex(gx)]
    o = isfinite.(v)
    count(o) >= 3 || continue
    c = get(cgrad(:viridis), (j - 1) / max(1, length(gy) - 1))
    lines!(ax2, gx[o], v[o]; color = c, linewidth = 2.2,
           label = @sprintf("sigma_gzz = %.2f", sg))
    scatter!(ax2, gx[o], v[o]; color = c, markersize = 7)
end
vlines!(ax2, [3.44]; color = (:orange, 0.9), linestyle = :dash, linewidth = 2)
vlines!(ax2, [3.8]; color = (:red, 0.6), linestyle = :dot, linewidth = 2)
text!(ax2, 3.45, 40; text = "published", color = :orange, fontsize = 10, rotation = pi/2)
text!(ax2, 3.81, 40; text = "by-eye", color = :red, fontsize = 10, rotation = pi/2)
scatter!(ax2, [gbest[2]], [gbest[1]]; marker = :star5, markersize = 22, color = :white,
         strokecolor = :black, strokewidth = 1.5)
axislegend(ax2; position = :rt, labelsize = 9, nbanks = 2)

# ---------------------------------------------------------------- 3  K/M surface
if !isempty(s2)
    jx, jy, JZ = surface(s2, "J1_meV", "sigma_J")
    ax3 = Axis(fig[1, 3], xlabel = "J1 (meV)", ylabel = "sigma_J",
               title = "K and M only: chi2 rises with J1, falls with sigma_J")
    hm3 = heatmap!(ax3, jx, jy, log10.(JZ); colormap = :magma)
    Colorbar(fig[1, 3][1, 2], hm3, label = "log10 chi2_red")
    jb = argmin_grid(jx, jy, JZ)
    scatter!(ax3, [jb[2]], [jb[3]]; marker = :star5, markersize = 26, color = :white,
             strokecolor = :black, strokewidth = 1.5)
    scatter!(ax3, [0.25], [0.5]; marker = :xcross, markersize = 18, color = :cyan,
             strokecolor = :black, strokewidth = 1)
    text!(ax3, 0.255, 0.53; text = "by-eye", color = :cyan, fontsize = 11)
    text!(ax3, jx[1] + 0.004, jy[end] - 0.06;
          text = "minimum on BOTH edges\n=> running away", color = :white, fontsize = 10)
    @printf("K/M minimum: J1 = %.3f, sigma_J = %.2f, chi2_red = %.1f\n", jb[2], jb[3], jb[1])

    # ------------------------------------------------------------ 4  sigma_J slices
    ax4 = Axis(fig[2, 1], xlabel = "sigma_J", ylabel = "chi2_red (K and M)",
               title = "sigma_J lowers chi2 at EVERY J1 -- it absorbs misfit")
    for (i, j1) in enumerate(jx)
        v = [JZ[i, k] for k in eachindex(jy)]
        o = isfinite.(v)
        c = get(cgrad(:plasma), (i - 1) / max(1, length(jx) - 1))
        lines!(ax4, jy[o], v[o]; color = c, linewidth = 2.5,
               label = @sprintf("J1 = %.2f", j1))
        scatter!(ax4, jy[o], v[o]; color = c, markersize = 8)
    end
    axislegend(ax4; position = :rt, labelsize = 10)
    text!(ax4, 0.02, 250; text = "no crossing: the preferred sigma_J\n" *
          "does not depend on J1", fontsize = 10, color = :black)
end

# ---------------------------------------------------------------- 5  stage 4 trajectory
if !isempty(s4)
    ax5 = Axis(fig[2, 2], xlabel = "stage-4 trial", ylabel = "chi2_red (all six cuts)",
               title = "Joint refinement: the factorized point is only a waypoint",
               yscale = log10)
    ch = [num(r, "chi2_red") for r in s4]
    acc = [lowercase(get(r, "accepted", "false")) == "true" for r in s4]
    par = [get(r, "parameter", "?") for r in s4]
    xs = 1:length(ch)
    lines!(ax5, xs, ch; color = (:grey, 0.6), linewidth = 1.5)
    scatter!(ax5, xs[acc], ch[acc]; color = :seagreen, markersize = 13,
             label = "accepted")
    scatter!(ax5, xs[.!acc], ch[.!acc]; color = :indianred, marker = :xcross,
             markersize = 11, label = "rejected")
    # Running best, which is what the search actually follows.
    best = accumulate(min, ch)
    lines!(ax5, xs, best; color = :black, linewidth = 2.5, linestyle = :dash,
           label = "running best")
    for (i, p) in enumerate(par)
        acc[i] || continue
        text!(ax5, i, ch[i] * 1.18; text = p, fontsize = 9, align = (:center, :bottom))
    end
    axislegend(ax5; position = :rt, labelsize = 10)
    @printf("stage 4: chi2_red %.4g -> %.4g over %d trials\n", ch[1], minimum(ch), length(ch))
end

# ---------------------------------------------------------------- 6  per-cut chi2
pc = filter(k -> startswith(k, "chi2_0"), s2hdr)
if !isempty(pc) && !isempty(s2)
    jx, jy, JZ = surface(s2, "J1_meV", "sigma_J")
    jb = argmin_grid(jx, jy, JZ)
    findrow(rows, j1, sj) = begin
        hit = nothing
        for r in rows
            if num(r, "J1_meV") ≈ j1 && num(r, "sigma_J") ≈ sj
                hit = r
                break
            end
        end
        hit
    end
    row = findrow(s2, jb[2], jb[3])
    ax6 = Axis(fig[2, 3], ylabel = "chi2_red", xticks = (1:length(pc) + 2,
               vcat(["(0,1,0)\n9 T", "(0,1,0)\n14 T"],
                    [replace(replace(k, "chi2_" => ""), "_" => " ") for k in pc])),
               title = "Zone centre fits ~25x better than the dispersive cuts
" *
                       "(Gamma from stage 1 at its own optimum; K/M from stage 2 at its own)",
               titlesize = 12,
               yscale = log10, xticklabelsize = 9)
    # Gamma values come from the stage-1 surface at its own minimum; K/M from stage 2.
    gvals = [gbest[1], gbest[1]]
    kvals = row === nothing ? fill(NaN, length(pc)) : [num(row, k) for k in pc]
    vals = vcat(gvals, kvals)
    cols = vcat(fill(:steelblue, 2), fill(:indianred, length(pc)))
    barplot!(ax6, 1:length(vals), vals; color = cols)
    for (i, v) in enumerate(vals)
        isfinite(v) || continue
        text!(ax6, i, v * 1.08; text = @sprintf("%.0f", v), fontsize = 10,
              align = (:center, :bottom))
    end
    text!(ax6, 1.0, maximum(finite(vals)) * 0.35;
          text = "blue: Zeeman sector only\n(gzz, sigma_gzz) -- model works\n\n" *
                 "red: exchange sector\n-- model fails here",
          fontsize = 10, align = (:left, :top))
end

out = joinpath(FDIR, "gamma_first_scan_diagnostics.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
