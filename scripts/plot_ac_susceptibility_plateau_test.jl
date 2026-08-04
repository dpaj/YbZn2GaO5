#!/usr/bin/env julia
# AC susceptibility at ~20 mK to 18 T: what the NHMFL SCM1 dataset can and cannot tell us.
# DATA ONLY -- no model, no Sunny, ~1 min.
#
#   julia --project=. scripts/plot_ac_susceptibility_plateau_test.jl
#
# WHAT THIS WAS BUILT TO TEST, AND WHY IT COULD NOT BE DONE AS INTENDED
#
# The intended test was sharp: the published model wants field-induced ordered phases with a
# magnetisation PLATEAU, a plateau is a DIP TOWARD ZERO in dM/dH, and AC susceptibility measures
# chi' = dM/dH directly at 20 mK to 18 T -- beyond both the 14 T DC data and the 14 T neutron
# measurements, cold enough that thermal rounding cannot hide a feature, and scale-free so the
# absence of absolute units costs nothing.
#
# Reading the run log (SCM1_July2025.xlsx) rather than guessing from the column names changes what
# is on offer. The probe carried THREE pickup coils and only one held our sample:
#
#   T3   YbZn2GaO5, "para" = B || c     <-- OURS, and the ONLY one
#   T1   LuCu(OH)Br                     <-- a DIFFERENT COMPOUND, another group on the same probe
#   B1   not listed in the log          <-- yet carries the LARGEST signal, 6x T3
#
# So: THERE IS NO PERPENDICULAR AC MEASUREMENT. An earlier version of this script assumed B1 was
# the perpendicular crystal, which is wrong -- and the assumption was self-refuting, because B1's
# in-phase channel crosses zero and goes NEGATIVE above ~12.4 T, which no sample susceptibility can
# do. Whatever B1 is, it is dominated by instrumental field dependence.
#
# ALL THREE COILS ARE READ ON LOCK-IN 1. The three SR830s sit at 991, 313 and 137 Hz, but only 991
# Hz drives in these runs, so the `x2`/`x3`/`y2`/`y3` columns are noise at 1e-9, three orders down.
# Use `x1`,`y1` only. Magnitude sqrt(x^2+y^2) is reported throughout because it is free of the
# lock-in phase convention, and the phase is plotted separately as a diagnostic.
#
# THE RESULT IS A CONTRADICTION WITH OUR OWN DC DATA, NOT A PLATEAU MEASUREMENT
#
# T3's magnitude rises to 1 T, is FLAT from ~1 to ~10 T, then falls, reaching ~40% at 18 T. The DC
# dM/dH falls to ~5% of its 1 T value by 10 T, because the sample saturates near 5 T. That is a
# factor ~8 disagreement at 10 T, in the wrong direction for temperature to explain: cooling from
# 2.5 K to 20 mK SHARPENS a saturation, it cannot flatten one.
#
# TWO BACKGROUND MODELS WERE TRIED AND BOTH FAIL.
#
#   A CONSTANT COMPLEX COIL OFFSET. Subtracting the full 18 T complex value leaves 6.77e-7 at 2 T
#   against 6.63e-7 at 8 T -- still flat. The flat region survives the obvious background model.
#
#   DIFFERENCING THE TWO TEMPERATURES. Runs 047/048 are at 450 mK (the log's "4.6 mA" is a heater
#   current), so 015-047 should cancel any term that depends on coil and magnet but not on the spin
#   state. The difference is instead at the noise level from 2-7 T and GROWS to 7.8e-8 V by 14-18 T,
#   largest exactly where the sample is most saturated, and it changes sign below 2 T. No sample
#   susceptibility does that, so the non-sample term is not temperature-independent either.
#
# So the AC magnitude CANNOT be read as chi'(sample) with what is in this dataset. What is missing is
# specific and experimental, not analytical: an empty-coil run at matching field and temperature,
# the coil constants, and an account of what coil B1 held. The DC data are the better characterised
# side of the contradiction by a wide margin -- absolute units, two instruments, two crystals, and
# agreement with the published SI to four digits -- so the contradiction is not evidence against them.
#
# WHAT SURVIVES ANYWAY, AND IT IS WORTH HAVING
#
# One class of statement is immune to any smooth background, because a smooth background cannot
# create or cancel a sharp feature: there is NO local dip, step, spike or hysteresis anywhere in
# 0-18 T, in either sweep direction, in two independent sweep pairs (015/016 and 047/048). A
# plateau EDGE or a first-order transition would show up as exactly such a feature. That is a
# weaker negative result than the one intended -- it bounds sharp features, not the overall
# magnitude of chi' -- but it is a real one, and it does not depend on the calibration that is
# missing.

using Printf, Statistics, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const DIR = SV.sv_repo_path(REPO_ROOT, "data/ac_susceptibility/nhmfl")
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/ac_plateau_test")
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/ac_plateau_test")
mkpath(FDIR); mkpath(TDIR)

# Coil -> what is actually in it, from SCM1_July2025.xlsx. Never infer this from column names.
const COILS = ["T3" => "YbZn2GaO5, B ∥ c  (OURS)",
               "T1" => "LuCu(OH)Br  (different compound)",
               "B1" => "unlisted in run log  (reference / instrumental)"]
# 047/048 are NOT a repeat of 015/016. The run log's "4.6 mA" is a heater current, and Tmc confirms
# it: 015/016 sit at 20 mK, 047/048 at 450-500 mK. That makes them a TEMPERATURE PAIR, which is
# what rescues the dataset -- see the differencing argument in the header.
const SWEEPS = ["015" => "0 → 18 T, 20 mK", "016" => "18 → 0 T, 20 mK",
                "047" => "0 → 18 T, 450 mK", "048" => "18 → 0 T, 450 mK"]
const COLD, WARM = "015", "047"

"Read one SCM1 sweep: tab-separated, header on the line after the sequence history."
function read_scm1(run::AbstractString)
    path = joinpath(DIR, "Friday_SCM1_July2025.$(run).txt")
    lines = readlines(path)
    i = findfirst(l -> occursin("End of Sequence History", l), lines)
    i === nothing && error("No sequence-history marker in $path")
    hdr = strip.(String.(split(lines[i+1], '\t')))
    cols = Dict{String,Vector{Float64}}(h => Float64[] for h in hdr if !isempty(h))
    for l in lines[i+2:end]
        f = String.(split(l, '\t'))
        length(f) < length(hdr) && continue
        vals = [tryparse(Float64, strip(x)) for x in f[1:length(hdr)]]
        any(isnothing, vals) && continue
        for (k, h) in enumerate(hdr)
            isempty(h) || push!(cols[h], vals[k])
        end
    end
    return cols
end

"Bin a field sweep onto a uniform grid; the sweep is step-and-hold so fields repeat."
function bin_sweep(H, V; edges = 0.0:0.1:18.0)
    c, m, s = Float64[], Float64[], Float64[]
    for k in 1:(length(edges)-1)
        sel = (H .>= edges[k]) .& (H .< edges[k+1])
        count(sel) == 0 && continue
        push!(c, 0.5*(edges[k]+edges[k+1])); push!(m, mean(V[sel]))
        push!(s, count(sel) > 2 ? std(V[sel]) : 0.0)
    end
    return c, m, s
end

"Magnitude and phase of one coil on lock-in 1, binned in field. Magnitude is phase-convention free."
function coil_trace(cols, run, coil)
    H = cols["Field_$(run)"]
    X = cols["$(coil)x1_$(run)"]; Y = cols["$(coil)y1_$(run)"]
    ok = isfinite.(H) .& isfinite.(X) .& isfinite.(Y) .& (H .>= -0.05) .& (H .<= 18.2)
    H, X, Y = H[ok], X[ok], Y[ok]
    mag = sqrt.(X.^2 .+ Y.^2)
    ph  = atand.(Y, X)
    cm, mm, ms = bin_sweep(H, mag)
    _,  pm, _  = bin_sweep(H, ph)
    return (; field = cm, mag = mm, mag_sd = ms, phase = pm)
end

data = Dict{Tuple{String,String},Any}()
for (run, _) in SWEEPS
    cols = read_scm1(run)
    T = cols["Tmc_$(run)"]
    @printf("run %s (%-18s): %5d rows, Tmc %.4f-%.4f K\n", run,
            Dict(SWEEPS)[run], length(cols["Field_$(run)"]), minimum(T), maximum(T))
    for (coil, _) in COILS
        data[(run, coil)] = coil_trace(cols, run, coil)
    end
end

"Normalise a trace to its own value nearest 1 T, so only shape is compared."
function norm1T(f, v)
    k = argmin(abs.(f .- 1.0))
    return v ./ (abs(v[k]) < eps() ? 1.0 : v[k])
end

# DC comparison: dM/dH from the 2.5 K DynaCool data, in absolute units.
function dc_derivative(rel, mass_mg)
    lines = readlines(SV.sv_repo_path(REPO_ROOT, rel))
    i = findfirst(l -> strip(l) == "[Data]", lines)
    hdr = strip.(String.(split(lines[i+1], ',')))
    ci = Dict(h => k for (k, h) in enumerate(hdr))
    H, M = Float64[], Float64[]
    for l in lines[i+2:end]
        f = String.(split(l, ','))
        length(f) < length(hdr) && continue
        h = tryparse(Float64, f[ci["Magnetic Field (Oe)"]])
        m = tryparse(Float64, f[ci["Moment (emu)"]])
        (h === nothing || m === nothing) && continue
        # Positive branch only, so the centred difference cannot straddle the loop turning points.
        h/1e4 < 0.05 && continue
        push!(H, h/1e4); push!(M, m / (5585.0 * mass_mg*1e-3/453.53))
    end
    p = sortperm(H); H, M = H[p], M[p]
    d = similar(M)
    for k in eachindex(M)
        lo = max(firstindex(M), k-8); hi = min(lastindex(M), k+8)
        d[k] = (M[hi]-M[lo]) / max(1e-9, H[hi]-H[lo])
    end
    return H, d
end
Hdc, Ddc = dc_derivative("data/magnetization/ppms_2p5K/YZGO_BparaC_4.81MG_2.5K_06242026.DAT", 4.81)
Ddc_n = Ddc ./ Ddc[argmin(abs.(Hdc .- 1.0))]

fig = Figure(size = (1560, 1080))
Label(fig[0, 1:2],
      "NHMFL SCM1 AC susceptibility, ~20 mK to 18 T. Only coil T3 holds YbZn2GaO5, and only B ∥ c.";
      fontsize = 16, font = :bold)

# --- 1. All three coils, correctly attributed -----------------------------------------
ax1 = Axis(fig[1, 1], xlabel = "field (T)", ylabel = "|x + iy| at 991 Hz  (V)",
           title = "All three pickup coils — attribution from the run log, not the column names")
for ((coil, lbl), col) in zip(COILS, (:crimson, :seagreen, :grey45))
    t = data[("015", coil)]
    lines!(ax1, t.field, t.mag; color = col, linewidth = 2.6, label = "$coil: $lbl")
end
axislegend(ax1; position = :rt, labelsize = 9); xlims!(ax1, 0, 18.2)

# --- 2. Phase: where is the measurement even stable? ----------------------------------
ax2 = Axis(fig[1, 2], xlabel = "field (T)", ylabel = "phase (deg)",
           title = "Lock-in phase — T3 is stable to ~4° out to 12 T; B1 rotates 79°")
for ((coil, _), col) in zip(COILS, (:crimson, :seagreen, :grey45))
    t = data[("015", coil)]
    lines!(ax2, t.field, t.phase; color = col, linewidth = 2.6, label = coil)
end
vlines!(ax2, [12.0]; color = (:black, 0.45), linestyle = :dash)
text!(ax2, 12.3, 60; text = "T3 phase begins to\nrotate above ~12 T", fontsize = 10, color = :grey25)
axislegend(ax2; position = :lt, labelsize = 9); xlims!(ax2, 0, 18.2)

# --- 3. The contradiction with our own DC data ----------------------------------------
ax3 = Axis(fig[2, 1], xlabel = "field (T)", ylabel = "normalised to own value at 1 T",
           title = "T3 against DC dM/dH — a factor ~8 disagreement at 10 T")
t = data[("015", "T3")]
lines!(ax3, t.field, norm1T(t.field, t.mag); color = :crimson, linewidth = 2.8,
       label = "AC |T3| , 20 mK, B ∥ c")
lines!(ax3, Hdc, Ddc_n; color = :seagreen, linewidth = 2.6, linestyle = :dash,
       label = "DC dM/dH, 2.5 K, B ∥ c")
vlines!(ax3, [10.0]; color = (:black, 0.35), linestyle = :dot)
hlines!(ax3, [0.0]; color = (:black, 0.35), linestyle = :dash)
text!(ax3, 4.2, 0.60;
      text = "Cooling 2.5 K → 20 mK SHARPENS a saturation.\n" *
             "It cannot flatten one, so temperature does not\n" *
             "reconcile these. Subtracting the 18 T complex\n" *
             "value as a coil offset leaves the flat region intact.",
      fontsize = 10, color = :grey20)
axislegend(ax3; position = :rt, labelsize = 9); xlims!(ax3, 0, 18.2)

# --- 4. AN ATTEMPTED RESCUE, AND WHY IT FAILS -----------------------------------------
# The idea: the non-sample contribution should be a property of the coil and the magnet, not of the
# sample's spin state, hence the same at 20 mK and 450 mK at any given field. Differencing would then
# cancel it with no coil constant and no empty-coil run needed, leaving the temperature-DEPENDENT
# part of the sample response -- which MUST collapse to zero once the sample is saturated.
#
# IT DOES NOT WORK, and the way it fails is the useful part. The difference is at the noise level
# from 2 to 7 T but then GROWS monotonically to 7.8e-8 V by 14-18 T, i.e. it is largest exactly
# where the sample is most completely saturated and a sample response must be smallest. It also
# changes sign below 2 T. No sample susceptibility behaves that way, so the non-sample term is NOT
# temperature-independent, and this panel is evidence for that conclusion rather than a measurement
# of chi'. (Tcoil differs between the runs, and 048 drifts to 0.60 K, so there are candidate
# mechanisms; none of them can be pinned down without the missing calibration data.)
cold, warm = data[(COLD, "T3")], data[(WARM, "T3")]
common = [(k, argmin(abs.(warm.field .- cold.field[k]))) for k in eachindex(cold.field)]
common = [(k, j) for (k, j) in common if abs(warm.field[j] - cold.field[k]) < 0.06]
Bd  = [cold.field[k] for (k, _) in common]
dif = [cold.mag[k] - warm.mag[j] for (k, j) in common]
sd  = [sqrt(cold.mag_sd[k]^2 + warm.mag_sd[j]^2) for (k, j) in common]

ax4 = Axis(fig[2, 2], xlabel = "field (T)", ylabel = "|T3|(20 mK) − |T3|(450 mK)   (V)",
           title = "Attempted rescue by differencing temperatures — it FAILS, informatively")
band!(ax4, Bd, dif .- sd, dif .+ sd; color = (:purple, 0.20))
lines!(ax4, Bd, dif; color = :purple, linewidth = 2.8, label = "AC temperature difference")
hlines!(ax4, [0.0]; color = (:black, 0.5), linestyle = :dash)
noise = median(sd[Bd .>= 8.0])
hspan!(ax4, -2noise, 2noise; color = (:grey, 0.18))
ax4b = Axis(fig[2, 2], yaxisposition = :right, ylabel = "DC dM/dH  (μB / T / Yb)",
            ylabelcolor = :seagreen, yticklabelcolor = :seagreen)
hidespines!(ax4b); hidexdecorations!(ax4b); linkxaxes!(ax4, ax4b)
lines!(ax4b, Hdc, Ddc; color = :seagreen, linewidth = 2.4, linestyle = :dash,
       label = "DC dM/dH, 2.5 K")
for a in (ax4, ax4b); xlims!(a, 0, 18.2); end
axislegend(ax4; position = :rt, labelsize = 9)
text!(ax4, 2.6, maximum(dif)*0.42;
      text = "If this were the sample, it would track the green DC\n" *
             "curve and die at saturation. Instead it is at the noise\n" *
             "level from 2–7 T and LARGEST at 14–18 T, where the\n" *
             "sample is most saturated. So the non-sample term is not\n" *
             "temperature-independent, and differencing cannot\n" *
             "substitute for the missing calibration.\n" *
             "Grey band = ±2× the scatter above 8 T.",
      fontsize = 9.5, color = :grey20)

Label(fig[3, 1:2],
      "WHAT SURVIVES: a smooth instrumental background cannot create or cancel a SHARP feature, so " *
      "the absence of any dip, step, spike or up/down hysteresis anywhere in 0–18 T bounds " *
      "first-order transitions and plateau edges even though the magnitude of chi' is uncalibrated. " *
      "WHAT DOES NOT: chi' itself, because both candidate background models fail (see panel 4). " *
      "TO FIX IT, three things are needed from the experiment, not from analysis — an empty-coil run " *
      "at matching field and temperature, the coil constants, and an account of what coil B1 held.";
      fontsize = 10, color = :grey30, word_wrap = true, tellwidth = false)

open(joinpath(TDIR, "ac_coil_magnitudes.csv"), "w") do io
    println(io, "run,sweep,coil,contents,field_T,magnitude_V,magnitude_sd,phase_deg")
    for (run, slbl) in SWEEPS, (coil, clbl) in COILS
        tt = data[(run, coil)]
        for k in eachindex(tt.field)
            @printf(io, "%s,%s,%s,%s,%.3f,%.6e,%.6e,%.3f\n", run, slbl, coil,
                    replace(clbl, "," => ";"), tt.field[k], tt.mag[k], tt.mag_sd[k], tt.phase[k])
        end
    end
end

# Quantify "no sharp feature" instead of eyeballing it: compare the largest second difference of
# the T3 shape against the sweep-to-sweep reproducibility over the same field range.
println("\nSharp-feature bound on T3 (YbZn2GaO5, B ∥ c):")
ref = data[("015", "T3")]; rn = norm1T(ref.field, ref.mag)
d2 = [abs(rn[k-1] - 2rn[k] + rn[k+1]) for k in 2:(length(rn)-1)]
sel = ref.field[2:end-1] .>= 1.0
pairs = [("015", "016"), ("047", "048")]
repro = Float64[]
for (a, b) in pairs
    ta, tb = data[(a, "T3")], data[(b, "T3")]
    na, nb = norm1T(ta.field, ta.mag), norm1T(tb.field, tb.mag)
    for k in eachindex(ta.field)
        j = argmin(abs.(tb.field .- ta.field[k]))
        abs(tb.field[j] - ta.field[k]) < 0.06 && ta.field[k] >= 1.0 && push!(repro, abs(na[k]-nb[j]))
    end
end
@printf("  largest curvature |d2(chi')| above 1 T   : %.4f per (0.1 T)^2\n", maximum(d2[sel]))
@printf("  at field                                 : %.2f T\n",
        ref.field[2:end-1][sel][argmax(d2[sel])])
@printf("  up/down reproducibility, median          : %.4f   (max %.4f)\n",
        median(repro), maximum(repro))
@printf("  monotonic decrease above 2 T?            : %s\n",
        all(diff(rn[ref.field .>= 2.0]) .<= 0.01) ? "yes, within 0.01" : "no")
for run in (COLD, WARM)
    tt = data[(run, "T3")]; nn = norm1T(tt.field, tt.mag)
    k = argmin(nn[tt.field .>= 2.0])
    @printf("  run %s: minimum of chi' above 2 T at %.2f T (%s)\n", run,
            tt.field[tt.field .>= 2.0][k],
            tt.field[tt.field .>= 2.0][k] > 17.5 ? "the ENDPOINT -- no interior dip" : "INTERIOR DIP")
end

println("\nTemperature-differenced sample response, |T3|(20 mK) - |T3|(450 mK):")
@printf("  peak                     : %.3e V at %.2f T\n", maximum(dif), Bd[argmax(dif)])
@printf("  scatter above 8 T (noise): %.3e V\n", noise)
for B in (1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 14.0, 18.0)
    k = argmin(abs.(Bd .- B))
    j = argmin(abs.(Hdc .- B))
    @printf("  %5.1f T : difference %+.3e V  (%6.1f%% of peak)   DC dM/dH %.4f uB/T (%5.1f%% of its 1 T value)\n",
            B, dif[k], 100*dif[k]/maximum(dif), Hdc[j] <= maximum(Hdc) ? Ddc[j] : NaN,
            100*Ddc_n[j])
end
kfall = findfirst(k -> Bd[k] > Bd[argmax(dif)] && abs(dif[k]) < 2noise, eachindex(Bd))
@printf("  falls into the noise band above : %s\n",
        kfall === nothing ? "never within 18 T" : @sprintf("%.2f T", Bd[kfall]))

out = joinpath(FDIR, "ac_susceptibility_plateau_test.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
