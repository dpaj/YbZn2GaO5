#!/usr/bin/env julia
# STAGE 0 of the background programme: characterise the Ei = 3.32 meV background from its own
# field set, before anything is built on top of it. DATA ONLY -- no model, no Sunny, ~1 min.
#
#   julia --project=. scripts/background_stage0_ei332.jl
#
# WHY THIS HAS TO COME FIRST
#
# The plan is to use the Ei = 3.32 data to estimate the true signal across the region where
# the Ei = 4.65 background construction has to interpolate. That only works if the 3.32 data
# have had THEIR OWN background removed first -- they carry extrinsic spectrometer effects of
# their own, notably a time-independent pile-up on the high-energy side near the kinematic
# limit. Subtracting raw 3.32 from raw 4.65 gives the DIFFERENCE of two backgrounds, not the
# 4.65 background.
#
# THE LOGIC, AND WHY IT IS TESTED RATHER THAN ASSUMED
#
# The magnet background sits at a FIXED energy while the signal MOVES WITH FIELD. So at each
# energy one takes whichever field has its signal furthest away -- which is what
# min-over-fields does. The low-energy side is therefore best estimated at 14 T, where the
# signal has moved up and vacated it.
#
# But "signal-free" cannot be assumed anywhere. The 0 T data are NOT signal-free: there is
# still a magnon at lower energy plus a diffuse zero-field continuum whose extent is not
# established, and modelling it with KPM is the hardest regime in this project -- at zero
# field the system is neither saturated nor near it, so the ground state is a complicated
# texture that may not be unique.
#
# So instead of declaring windows, this measures a statistic that says where the field-sweep
# logic actually works: the background is FIELD-INDEPENDENT, so wherever two or more fields
# AGREE, that common value is plausibly background. Where only the single lowest field is low,
# the others still carry signal and the minimum is an upper bound, not an estimate. The gap
# between the lowest and second-lowest value across fields is therefore a direct, per-energy
# measure of how much to trust the background there -- and it is exactly the quantity that
# should later become a variance.

using Printf, Statistics, LinearAlgebra, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const controls = SV.sv_load_controls(REPO_ROOT)
const DIR = SV.sv_repo_path(REPO_ROOT, controls["paths"]["neutron_1d_dir"])
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/background_stage0")
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/background_stage0")
mkpath(FDIR); mkpath(TDIR)

const EI = Float64(get(ENV, "YZGO_STAGE0_EI", "") == "" ? 3.32 :
                   parse(Float64, ENV["YZGO_STAGE0_EI"]))
const TBASE = 0.07

scans = Dict{Tuple{Float64,String},Any}()
for f in sort(readdir(DIR))
    endswith(lowercase(f), ".dat") || continue
    path = joinpath(DIR, f)
    meta = try SV.sv_parse_neutron_1d_filename(path) catch; continue end
    (meta.Ei_meV ≈ EI && meta.temperature_K ≈ TBASE) || continue
    scans[(meta.field_T, meta.qtag)] = SV.sv_load_neutron_raw_scan_1d(path, controls)
end
isempty(scans) && error("No Ei = $EI, T = $TBASE scans found in $DIR")

qtags = sort(unique(String[k[2] for k in keys(scans)]))
fields = sort(unique(Float64[k[1] for k in keys(scans)]))
@printf("Ei = %.2f meV, T = %.2f K: %d qtags %s, %d fields %s\n",
        EI, TBASE, length(qtags), join(qtags, ","), length(fields), string(fields))

qlabel(t) = (p = split(t, "_"); length(p) == 3 ?
             "(" * join(replace.(p, "p" => "."), ", ") * ")" : t)

"""
Per-energy order statistics across fields. `lo` is the min-over-fields background estimate,
`lo2` the second-lowest. `agree = (lo2 - lo) / lo` is small where at least two fields give
the same value, which is the condition under which the minimum is a background ESTIMATE
rather than merely an upper bound.
"""
function field_stats(byfield, fields, E)
    n = length(E)
    lo = fill(NaN, n); lo2 = fill(NaN, n); hi = fill(NaN, n); which = fill(NaN, n)
    for i in 1:n
        vals = Float64[]; fs = Float64[]
        for B in fields
            s = byfield[B]
            y = SV.sv_interp1(s.energy_meV, s.intensity, [E[i]])[1]
            isfinite(y) || continue
            push!(vals, y); push!(fs, B)
        end
        length(vals) >= 2 || continue
        p = sortperm(vals)
        lo[i] = vals[p[1]]; lo2[i] = vals[p[2]]; hi[i] = vals[p[end]]
        which[i] = fs[p[1]]
    end
    agree = (lo2 .- lo) ./ max.(abs.(lo), eps())
    return (; lo, lo2, hi, which, agree)
end

rows = Vector{NTuple{8,Any}}()
fig = Figure(size = (620 * length(qtags) + 60, 1080))
Label(fig[0, 1:length(qtags)],
      "Stage 0 -- Ei = $(EI) meV background from its own field set.  " *
      "Background is field-INDEPENDENT, so agreement between fields is the test.";
      fontsize = 16, font = :bold)

for (col, q) in enumerate(qtags)
    byfield = Dict(B => scans[(B, q)] for B in fields if haskey(scans, (B, q)))
    length(byfield) >= 2 || continue
    fs = sort(collect(keys(byfield)))
    E = byfield[fs[1]].energy_meV
    st = field_stats(byfield, fs, E)

    # ---- row 1: the three fields, and the min-over-fields estimate ----------------
    ax1 = Axis(fig[1, col], ylabel = col == 1 ? "raw intensity (arb.)" : "",
               title = "$(qlabel(q))   fields overlaid")
    cols = [:steelblue, :seagreen, :indianred, :purple]
    for (k, B) in enumerate(fs)
        s = byfield[B]
        lines!(ax1, s.energy_meV, s.intensity; color = cols[mod1(k, length(cols))],
               linewidth = 2.0, label = @sprintf("%.0f T", B))
        scatter!(ax1, s.energy_meV, s.intensity; color = cols[mod1(k, length(cols))],
                 markersize = 5)
    end
    lines!(ax1, E, st.lo; color = :black, linewidth = 3.0, linestyle = :dash,
           label = "min over fields")
    ylims!(ax1, -0.0003, 1.15 * maximum(filter(isfinite, st.hi[E .>= 0.35])))
    xlims!(ax1, 0.0, 3.2)
    col == 1 && axislegend(ax1; position = :rt, labelsize = 10)

    # ---- row 2: do the fields AGREE? ---------------------------------------------
    ax2 = Axis(fig[2, col], ylabel = col == 1 ? "(2nd lowest - lowest) / lowest" : "",
               title = "agreement between fields (low = min is a real estimate)",
               yscale = log10)
    a = copy(st.agree)
    a[.!isfinite.(a) .| (a .<= 0)] .= NaN
    lines!(ax2, E, a; color = :black, linewidth = 2.0)
    scatter!(ax2, E, a; color = :black, markersize = 5)
    hlines!(ax2, [0.10]; color = :seagreen, linestyle = :dash, linewidth = 2)
    hlines!(ax2, [0.50]; color = :indianred, linestyle = :dash, linewidth = 2)
    text!(ax2, 0.05, 0.105; text = "10% -- two fields agree, trust the minimum",
          color = :seagreen, fontsize = 9)
    text!(ax2, 0.05, 0.52; text = "50% -- only one field is low, minimum is an upper bound",
          color = :indianred, fontsize = 9)
    xlims!(ax2, 0.0, 3.2)

    # ---- row 3: which field supplies the minimum, and the high-E pile-up ---------
    ax3 = Axis(fig[3, col], xlabel = "energy transfer (meV)",
               ylabel = col == 1 ? "field giving the minimum (T)" : "",
               title = "which field vacates each energy")
    scatter!(ax3, E, st.which; color = :darkorange, markersize = 7)
    ylims!(ax3, -1.5, maximum(fs) + 1.5)
    xlims!(ax3, 0.0, 3.2)

    # ---- report ------------------------------------------------------------------
    good = isfinite.(st.agree) .& (st.agree .< 0.10) .& (E .>= 0.3)
    @printf("\n%s (%s):\n", q, qlabel(q))
    @printf("  energies where >=2 fields agree within 10%%: %d of %d above 0.3 meV\n",
            count(good), count(E .>= 0.3))
    if any(good)
        eg = E[good]
        # Report contiguous runs, which are the usable windows.
        runs = Vector{Tuple{Float64,Float64}}(); s0 = eg[1]; prev = eg[1]
        for e in eg[2:end]
            if e - prev > 0.12
                push!(runs, (s0, prev)); s0 = e
            end
            prev = e
        end
        push!(runs, (s0, prev))
        @printf("  usable windows: %s\n",
                join([@sprintf("[%.2f, %.2f]", r[1], r[2]) for r in runs], "  "))
    end
    # High-energy pile-up: is the rise near the kinematic limit the SAME at all fields?
    # If it is, it is instrumental and subtractable; if not, it is not purely instrumental.
    hiE = E .>= 2.5
    if count(hiE) > 3
        sp = filter(isfinite, st.agree[hiE])
        @printf("  high-E (>2.5 meV) median field disagreement: %.1f%%  => %s\n",
                100 * median(sp),
                median(sp) < 0.15 ? "consistent across fields, i.e. INSTRUMENTAL" :
                                    "field-dependent, so NOT purely instrumental")
    end
    for i in eachindex(E)
        push!(rows, (q, E[i], st.lo[i], st.lo2[i], st.hi[i], st.which[i], st.agree[i],
                     isfinite(st.agree[i]) && st.agree[i] < 0.10))
    end
end

open(joinpath(TDIR, "ei$(replace(string(EI), "." => "p"))_field_stats.csv"), "w") do io
    println(io, "qtag,energy_meV,min_over_fields,second_min,max_over_fields," *
                "field_of_min_T,disagreement,trusted")
    for r in rows
        @printf(io, "%s,%.6g,%.6g,%.6g,%.6g,%.1f,%.6g,%s\n", r...)
    end
end

Label(fig[4, 1:length(qtags)],
      "Row 2 is the key panel. The background is field-independent, so where two or more " *
      "fields agree the min-over-fields value is a background ESTIMATE; where only the " *
      "single lowest field is low, the others still carry signal and the minimum is only an " *
      "UPPER BOUND. That per-energy disagreement is the quantity that should become a " *
      "background variance downstream. Note the 0 T scans are NOT signal-free -- a magnon " *
      "persists at lower energy plus a diffuse zero-field continuum -- so no window is " *
      "assumed here, only measured.";
      fontsize = 10, color = :grey30)

out = joinpath(FDIR, "background_stage0_ei$(replace(string(EI), "." => "p")).png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
println("wrote " * joinpath(TDIR, "ei$(replace(string(EI), "." => "p"))_field_stats.csv"))
