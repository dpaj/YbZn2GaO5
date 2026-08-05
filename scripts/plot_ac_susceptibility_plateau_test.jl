#!/usr/bin/env julia
# AC susceptibility at 20 mK to 18 T, both orientations. DATA ONLY -- no model, no Sunny, ~1 min.
#
#   julia --project=. scripts/plot_ac_susceptibility_plateau_test.jl
#
# WHY THIS IS THE SHARPEST AVAILABLE TEST OF THE CENTRAL CLAIM
#
# The published model wants field-induced ordered phases with a magnetisation PLATEAU. A plateau is
# a DIP TOWARD ZERO in dM/dH, and AC susceptibility measures chi' = dM/dH DIRECTLY -- M is the
# derived quantity here, not the measured one. At 20 mK and up to 18 T this reaches past both the
# 14 T DC data and the 14 T neutron measurements, is cold enough that thermal rounding cannot hide a
# feature, and being a scale-free shape feature, costs nothing for having no absolute units.
#
# THE COIL KEY IS IN THE SPREADSHEET, COLUMNS E/F/G ROWS 1-2. Column names are COIL POSITIONS:
#
#   B1   YbZn2GaO5, "perp" = B PERPENDICULAR to c
#   T1   LuCu(OH)Br            <-- a DIFFERENT COMPOUND, another group sharing the probe
#   T3   YbZn2GaO5, "para" = B PARALLEL to c
#
# BOTH ORIENTATIONS OF OUR SAMPLE ARE PRESENT. An earlier version of this script claimed there was
# no perpendicular measurement. That was wrong and the error was OURS: a self-closing empty cell in
# the xlsx made our parser swallow the neighbouring cell, so B1's entry was read as a bare integer
# and mistaken for a row label. The supporting argument -- that B1 goes negative above ~12.4 T and
# so cannot be a sample -- was ALSO wrong: a raw lock-in quadrature contains the coil's own
# mutual-inductance background, which is large here and can dominate, so either quadrature may be
# negative.
#
# HENCE THE METHOD: WORK IN THE COMPLEX PLANE. The coil background and the sample response have
# different phases, so a single quadrature mixes them and even the bare magnitude can be
# non-monotonic while the underlying sample response is not. Use z = x1 + i*y1 throughout.
# Only the 991 Hz drive produces signal -- x2/x3/y2/y3 sit at 1e-9, three orders down.
#
# THE INSTRUMENTAL DRIFT, AND THE IN-SITU REFERENCE THAT BEATS IT
#
# Referencing each channel to its own 18 T value does NOT isolate the sample: all three coils,
# INCLUDING the different compound, then show a common near-linear ramp to zero, which no set of
# three different samples would produce. That common part is instrumental.
#
# T1 is the handle. LuCu(OH)Br sits on the same probe, in the same magnet, at the same temperature,
# and does not saturate over this range, so its complex field dependence is a template for the
# drift. Subtracting a fitted multiple of it from each YZGO channel gives, normalised at 1 T:
#
#   B (T)     perp AC 20 mK   perp DC 2.5 K   para AC 20 mK   para DC 2.5 K
#     2           0.913           0.922           0.844           0.882
#     3           0.817           0.820           0.772           0.752
#     5           0.345           0.691           0.724           0.452
#     8           0.109           0.371           0.590           0.125
#
# THE PERPENDICULAR CHANNEL WORKS: ~1% agreement with DC at 1-3 T, then saturating FASTER, which is
# the correct direction for 20 mK against a 125x higher temperature -- cooling sharpens a saturation.
# THE PARALLEL CHANNEL DOES NOT: it stays high where DC says the sample is saturated. Unexplained.
#
# CAVEATS, because the T1 scale is FITTED. It is fitted on 10-18 T assuming the sample is saturated
# there, so the high-field end of the perp agreement is partly circular; the 1-8 T comparison is not.
# T1 is a different sample in a different coil, so it is a template for the drift's SHAPE, not a
# calibrated background. The fitted scales are large and of opposite sign for the two channels, which
# is a further reason to treat this as provisional. An empty-coil run at matching field and
# temperature, plus the coil constants, is what would settle it.
#
# WHAT NEEDS NO CALIBRATION AT ALL: a smooth background cannot create or cancel a SHARP feature, so
# the absence of any dip, step, spike or up/down hysteresis anywhere in 0-18 T bounds first-order
# transitions and plateau edges regardless of all of the above.

using Printf, Statistics, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const DIR = SV.sv_repo_path(REPO_ROOT, "data/ac_susceptibility/nhmfl")
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/ac_plateau_test")
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/ac_plateau_test")
mkpath(FDIR); mkpath(TDIR)

# Coil -> what is actually in it, from SCM1_July2025.xlsx columns E/F/G rows 1-2.
# NEVER infer this from the column names; getting it wrong once already cost a wrong conclusion.
const COILS = ["B1" => "YbZn2GaO5, B ⟂ c  (OURS, perp)",
               "T1" => "LuCu(OH)Br  (different compound — instrumental reference)",
               "T3" => "YbZn2GaO5, B ∥ c  (OURS, para)"]
const REFCOIL = "T1"          # the non-saturating in-situ template for the instrumental drift
# 047/048 are NOT a repeat of 015/016. The run log's "4.6 mA" is a heater current, and Tmc confirms
# it: 015/016 sit at 20 mK, 047/048 at 450-500 mK -- so they are a TEMPERATURE pair.
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

"""
Binned COMPLEX response of one coil on lock-in 1. Complex, not magnitude: the coil background and
the sample response have different phases, so a single quadrature mixes them and even |z| can be
non-monotonic while the sample response is not.
"""
function coil_trace(cols, run, coil)
    H = cols["Field_$(run)"]
    X = cols["$(coil)x1_$(run)"]; Y = cols["$(coil)y1_$(run)"]
    ok = isfinite.(H) .& isfinite.(X) .& isfinite.(Y) .& (H .>= -0.05) .& (H .<= 18.2)
    H, X, Y = H[ok], X[ok], Y[ok]
    f, xr, xs = bin_sweep(H, X)
    _, yr, ys = bin_sweep(H, Y)
    z = complex.(xr, yr)
    return (; field = f, z, mag = abs.(z), phase = rad2deg.(angle.(z)),
              mag_sd = sqrt.(xs.^2 .+ ys.^2))
end

data = Dict{Tuple{String,String},Any}()
for (run, _) in SWEEPS
    cols = read_scm1(run)
    T = cols["Tmc_$(run)"]
    @printf("run %s (%-20s): %5d rows, Tmc %.4f-%.4f K\n", run,
            Dict(SWEEPS)[run], length(cols["Field_$(run)"]), minimum(T), maximum(T))
    for (coil, _) in COILS
        data[(run, coil)] = coil_trace(cols, run, coil)
    end
end

at(t, B) = argmin(abs.(t.field .- B))

"""
Sample term for `coil`, using `REFCOIL` as an in-situ template for the instrumental drift:

    s(B) = [z(B) - z(Bmax)] - alpha * [zref(B) - zref(Bmax)]

`alpha` is least squares over `fitrange`, where the sample is assumed saturated so the residual
should vanish. THAT ASSUMPTION IS THE WEAK POINT, and it makes the high-field end partly circular;
the low-field comparison does not depend on it. Referencing to `Bmax` alone is NOT enough -- all
three coils then show a common near-linear ramp to zero, the different compound included, so that
part is instrumental and has to be removed by the template.
"""
function sample_term(coil, run; Bmax = 18.0, fitrange = (10.0, 18.0))
    t, r = data[(run, coil)], data[(run, REFCOIL)]
    z0, r0 = t.z[at(t, Bmax)], r.z[at(r, Bmax)]
    num = den = 0.0
    for k in eachindex(t.field)
        (fitrange[1] <= t.field[k] <= fitrange[2]) || continue
        d = t.z[k] - z0; rr = r.z[at(r, t.field[k])] - r0
        num += real(d) * real(rr) + imag(d) * imag(rr)
        den += abs2(rr)
    end
    alpha = den > 0 ? num / den : 0.0
    s = [(t.z[k] - z0) - alpha * (r.z[at(r, t.field[k])] - r0) for k in eachindex(t.field)]
    return (; field = t.field, amp = abs.(s), alpha)
end

"Normalise to the value nearest 1 T, so only shape is compared."
norm1T(f, v) = (k = argmin(abs.(f .- 1.0)); v ./ (abs(v[k]) < eps() ? 1.0 : v[k]))

# DC dM/dH in absolute units, for both orientations.
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
        h/1e4 < 0.05 && continue      # positive branch only, so the difference cannot straddle it
        push!(H, h/1e4); push!(M, m / (5585.0 * mass_mg*1e-3/453.53))
    end
    p = sortperm(H); H, M = H[p], M[p]
    d = similar(M)
    for k in eachindex(M)
        lo = max(firstindex(M), k-8); hi = min(lastindex(M), k+8)
        d[k] = (M[hi]-M[lo]) / max(1e-9, H[hi]-H[lo])
    end
    # The centred difference goes one-sided at the ends, so the last ~1 T is unreliable; trim it.
    keep = H .<= (maximum(H) - 1.0)
    return H[keep], d[keep]
end
Hperp, Dperp = dc_derivative("data/magnetization/ppms_2p5K/YZGO_BperpC_7.7MG_2.5K_06242026.DAT", 7.7)
Hpara, Dpara = dc_derivative("data/magnetization/ppms_2p5K/YZGO_BparaC_4.81MG_2.5K_06242026.DAT", 4.81)

fig = Figure(size = (1560, 1080))
Label(fig[0, 1:2],
      "NHMFL SCM1 AC susceptibility, 20 mK to 18 T. Both orientations of YbZn2GaO5 are present; " *
      "coil T1 holds a different compound and serves as an in-situ instrumental reference.";
      fontsize = 15, font = :bold)

ax1 = Axis(fig[1, 1], xlabel = "field (T)", ylabel = "|x + iy| at 991 Hz  (V)",
           title = "All three coils, raw magnitude — attribution from the xlsx key")
for ((coil, lbl), c) in zip(COILS, (:dodgerblue, :seagreen, :crimson))
    t = data[("015", coil)]
    lines!(ax1, t.field, t.mag; color = c, linewidth = 2.6, label = "$coil: $lbl")
end
axislegend(ax1; position = :rt, labelsize = 9); xlims!(ax1, 0, 18.2)

ax2 = Axis(fig[1, 2], xlabel = "Re z  (V)", ylabel = "Im z  (V)",
           title = "Complex trajectory, 0 → 18 T — the background is a large complex offset")
for ((coil, _), c) in zip(COILS, (:dodgerblue, :seagreen, :crimson))
    t = data[("015", coil)]
    lines!(ax2, real.(t.z), imag.(t.z); color = c, linewidth = 2.2, label = coil)
    scatter!(ax2, [real(t.z[at(t, 0.0)])], [imag(t.z[at(t, 0.0)])]; color = c, markersize = 11)
    scatter!(ax2, [real(t.z[at(t, 18.0)])], [imag(t.z[at(t, 18.0)])]; color = c,
             marker = :xcross, markersize = 13)
end
axislegend(ax2; position = :lt, labelsize = 9)
text!(ax2, 0.03, 0.04; space = :relative, fontsize = 9.5, color = :grey25,
      text = "circle = 0 T, cross = 18 T. A single quadrature cuts\n" *
             "across these paths, which is why x1 alone can go\n" *
             "negative without that meaning the coil holds no sample.")

for (col, (coil, name, Hdc, Ddc)) in enumerate([("B1", "PERPENDICULAR, B ⟂ c", Hperp, Dperp),
                                                ("T3", "PARALLEL, B ∥ c", Hpara, Dpara)])
    st = sample_term(coil, COLD)
    ax = Axis(fig[2, col], xlabel = "field (T)",
              ylabel = col == 1 ? "normalised to own value at 1 T" : "",
              title = "$name  ($coil)")
    lines!(ax, st.field, norm1T(st.field, st.amp); color = :crimson, linewidth = 2.8,
           label = "AC 20 mK, T1-referenced")
    lines!(ax, Hdc, Ddc ./ Ddc[argmin(abs.(Hdc .- 1.0))]; color = :seagreen, linewidth = 2.6,
           linestyle = :dash, label = "DC dM/dH, 2.5 K")
    hlines!(ax, [0.0]; color = (:black, 0.35), linestyle = :dash)
    xlims!(ax, 0, 18.2); ylims!(ax, -0.1, 1.35)
    axislegend(ax; position = :rt, labelsize = 9)
    verdict = col == 1 ?
        "AGREES with DC to ~1% at 1–3 T, then saturates\nFASTER — the correct direction for 20 mK\n" *
        "against 2.5 K, since cooling sharpens a\nsaturation. Usable over 1–8 T." :
        "STAYS HIGH where DC says the sample is\nsaturated. Unexplained, and not fixed by the\n" *
        "T1 reference. Do not read chi' off this channel\nuntil an empty-coil run exists."
    text!(ax, 0.03, 0.28; space = :relative, fontsize = 9.5,
          color = col == 1 ? :grey20 : :firebrick, text = verdict)
    text!(ax, 0.03, 0.04; space = :relative, fontsize = 9, color = :grey45,
          text = "fitted T1 scale alpha = $(round(st.alpha; digits = 2))")
end

Label(fig[3, 1:2],
      "The T1 scale is FITTED on 10–18 T assuming the sample is saturated there, so the high-field " *
      "end is partly circular; the 1–8 T comparison is not. T1 is a different sample in a different " *
      "coil, hence a template for the drift's SHAPE, not a calibrated background. INDEPENDENT OF ALL " *
      "OF THAT: a smooth background cannot create or cancel a sharp feature, so the absence of any " *
      "dip, step, spike or up/down hysteresis in 0–18 T bounds first-order transitions and plateau " *
      "edges. An empty-coil run at matching field and temperature, plus the coil constants, is what " *
      "would make chi' absolute and settle the parallel channel.";
      fontsize = 10, color = :grey30, word_wrap = true, tellwidth = false)

open(joinpath(TDIR, "ac_sample_term.csv"), "w") do io
    println(io, "coil,contents,run,field_T,raw_mag_V,raw_phase_deg,sample_amp_V,sample_norm_1T")
    for (coil, clbl) in COILS, (run, _) in SWEEPS
        t = data[(run, coil)]
        st = coil == REFCOIL ? nothing : sample_term(coil, run)
        nn = st === nothing ? nothing : norm1T(st.field, st.amp)
        for k in eachindex(t.field)
            @printf(io, "%s,%s,%s,%.3f,%.6e,%.3f,%s,%s\n", coil, replace(clbl, "," => ";"), run,
                    t.field[k], t.mag[k], t.phase[k],
                    st === nothing ? "" : @sprintf("%.6e", st.amp[k]),
                    nn === nothing ? "" : @sprintf("%.6f", nn[k]))
        end
    end
end

println("\nT1-referenced sample term against DC dM/dH, each normalised at 1 T:")
@printf("  %5s | %14s %14s | %14s %14s\n",
        "B(T)", "perp AC 20mK", "perp DC 2.5K", "para AC 20mK", "para DC 2.5K")
sp, spa = sample_term("B1", COLD), sample_term("T3", COLD)
np_, npa = norm1T(sp.field, sp.amp), norm1T(spa.field, spa.amp)
dp = Dperp ./ Dperp[argmin(abs.(Hperp .- 1.0))]
da = Dpara ./ Dpara[argmin(abs.(Hpara .- 1.0))]
for Bt in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0)
    kp = argmin(abs.(sp.field .- Bt)); ka = argmin(abs.(spa.field .- Bt))
    jp = argmin(abs.(Hperp .- Bt));    ja = argmin(abs.(Hpara .- Bt))
    @printf("  %5.1f | %14.3f %14s | %14.3f %14s\n", Bt, np_[kp],
            abs(Hperp[jp] - Bt) < 0.4 ? @sprintf("%.3f", dp[jp]) : "-",
            npa[ka], abs(Hpara[ja] - Bt) < 0.4 ? @sprintf("%.3f", da[ja]) : "-")
end
@printf("\n  fitted T1 scale: alpha(perp) = %.3f, alpha(para) = %.3f\n", sp.alpha, spa.alpha)

# SHARP-FEATURE BOUND. This needs no calibration, but it has to be a test of SHARPNESS, not of where
# the minimum sits: B1's raw magnitude has a broad shallow minimum near 12 T that is simply the
# instrumental background turning over, and a "where is the minimum" test reports that as a dip. A
# plateau edge or a first-order transition is instead a feature LOCALISED in field, so the right
# statistic is the local curvature measured against the up/down sweep reproducibility -- both are
# insensitive to any smooth background, however large.
println("\nSharp-feature bound (calibration-free): local curvature vs up/down reproducibility.")
println("  A plateau edge or first-order transition is LOCALISED in field, so it must show up as")
println("  curvature well above the sweep-to-sweep scatter. Units: fraction of the 1 T value.")
for (coil, _) in COILS
    coil == REFCOIL && continue
    up, dn = data[(COLD, coil)], data[("016", coil)]
    nu = norm1T(up.field, up.mag)
    nd = norm1T(dn.field, dn.mag)
    # Reproducibility: |up - down| at matched fields, above 1 T.
    repro = Float64[]
    for k in eachindex(up.field)
        up.field[k] >= 1.0 || continue
        j = argmin(abs.(dn.field .- up.field[k]))
        abs(dn.field[j] - up.field[k]) < 0.06 && push!(repro, abs(nu[k] - nd[j]))
    end
    # Curvature: second difference over the 0.1 T grid, above 1 T.
    d2 = [abs(nu[k-1] - 2nu[k] + nu[k+1]) for k in 2:(length(nu)-1)]
    Bs2 = up.field[2:end-1]
    sel = Bs2 .>= 1.0
    kmax = argmax(d2[sel])
    noise = median(repro)
    @printf("  %-3s: max curvature %.4f at %5.2f T;  up/down median %.4f, max %.4f;  ratio %.1fx  %s\n",
            coil, d2[sel][kmax], Bs2[sel][kmax], noise, maximum(repro),
            d2[sel][kmax] / max(noise, eps()),
            d2[sel][kmax] < 3 * maximum(repro) ?
                "-- NO sharp feature above the noise" : "-- CANDIDATE, investigate")
end

out = joinpath(FDIR, "ac_susceptibility_plateau_test.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
