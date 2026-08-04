#!/usr/bin/env julia
# The MPMS3 centring correction, explained. DATA ONLY, no model, ~1 min.
#
#   julia --project=. scripts/plot_mpms3_centering_correction.jl
#
# THE PROBLEM THIS SOLVES
#
# The 0.42 K M(H) curve used throughout this repo was digitised from a figure and saturated near
# 1.12 uB/Yb, well below g_par*S = 1.72 at g_par = 3.44. It also sat BELOW the 2.5 K measurement at
# every field, which is thermodynamically impossible for a paramagnet -- colder must give more
# magnetisation. Why the scale was off had been an open question.
#
# THE CAUSE
#
# It is a data-reduction choice, not a measurement failure. The MPMS3 file reports the moment in two
# columns from two different fits to the SQUID response:
#
#   DC Moment Fixed Ctr  -- sample centre held at the nominal position
#   DC Moment Free Ctr   -- sample centre floated as a fit parameter
#
# The `Moment (emu)` column is EMPTY in this file, so whoever reduced it had to pick one, and the
# digitised curve came from Fixed Centring. The sample actually sat ~3 mm off the nominal centre, so
# the fixed fit was solving at the wrong position and under-read by a uniform factor of ~1.49.
#
# The file contains its own verdict on which column to trust:
#   DC Fixed Fit quality ~0.28  (poor)      DC Free Fit quality ~0.96  (good)
#   nominal centre 39.686 mm    vs          fitted centre ~36.65 mm
#
# THE CONFIRMATION
#
# Free centring gives 1.651 uB/Yb at 6.975 T. An entirely independent measurement on a different
# crystal and a different instrument -- DynaCool VSM at 2.5 K -- gives 1.652 at the same field, and
# the Bag et al. PRL Supplement reads ~1.65 at 7 T. Three independent numbers agreeing to four
# digits.
#
# ONE MORE TRAP IN THE SAME FILE, worth knowing for any future reduction: the `Temperature (K)`
# column reads 1.56 K, which is the MPMS3 CHAMBER. The sample temperature is in a column called
# `He3 temp` and reads 0.42-0.50 K, consistent with the filename and with the 0.42 K this repo's
# M(H) protocol has always assumed.

using Printf, Statistics, DelimitedFiles, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const MOL = 453.53      # g/mol, YbZn2GaO5
const NA_MUB = 5585.0   # emu/mol per mu_B per formula unit, = N_A * mu_B
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/mpms3_centering")
mkpath(FDIR)

"Read a Quantum Design .dat file: returns a Dict of column name => Vector{Float64}."
function read_qd(path)
    lines = readlines(path)
    i = findfirst(l -> strip(l) == "[Data]", lines)
    i === nothing && error("No [Data] section in $path")
    hdr = strip.(String.(split(lines[i+1], ',')))
    cols = Dict{String,Vector{Float64}}(h => Float64[] for h in hdr if !isempty(h))
    for l in lines[i+2:end]
        f = String.(split(l, ','))
        length(f) < length(hdr) && continue
        for (k, h) in enumerate(hdr)
            isempty(h) && continue
            v = tryparse(Float64, strip(f[k]))
            push!(cols[h], v === nothing ? NaN : v)
        end
    end
    return cols
end

emu_to_muB(mass_mg) = 1.0 / (NA_MUB * mass_mg * 1e-3 / MOL)

# ---- the 0.42 K MPMS3 file, both centring choices ---------------------------------
f04 = only(filter(p -> endswith(lowercase(p), ".dat"),
                  readdir(joinpath(REPO_ROOT, "data/magnetization/mpms3_0p4K"); join=true)))
c04 = read_qd(f04)
conv04 = emu_to_muB(12.2)
H04 = c04["Magnetic Field (Oe)"] ./ 1e4
Mfix = c04["DC Moment Fixed Ctr (emu)"] .* conv04
Mfree = c04["DC Moment Free Ctr (emu)"] .* conv04
Tsample = c04["He3 temp"]
Tchamber = c04["Temperature (K)"]
qfix, qfree = c04["DC Fixed Fit"], c04["DC Free Fit"]
ctr_fit, ctr_nom = c04["DC Calculated Center (mm)"], c04["Center Position (mm)"]

# ---- the independent 2.5 K DynaCool VSM measurement, different crystal -------------
f25 = joinpath(REPO_ROOT, "data/magnetization/ppms_2p5K/YZGO_BparaC_4.81MG_2.5K_06242026.DAT")
c25 = read_qd(f25)
H25 = c25["Magnetic Field (Oe)"] ./ 1e4
M25 = c25["Moment (emu)"] .* emu_to_muB(4.81)

# ---- the digitised curve that has been used until now ------------------------------
dig = readdlm(joinpath(REPO_ROOT, "data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv"),
               ',', String)
r0 = any(isnothing, tryparse.(Float64, dig[1, :])) ? 2 : 1
Hdig = [something(tryparse(Float64, s), NaN) for s in dig[r0:end, 1]]
Mdig = [something(tryparse(Float64, s), NaN) for s in dig[r0:end, 2]]

@printf("0.42 K MPMS3: sample T (He3 temp) = %.3f-%.3f K, chamber T = %.2f K\n",
        minimum(Tsample), maximum(Tsample), mean(Tchamber))
@printf("  fit quality: fixed %.3f, free %.3f     centre: nominal %.3f mm, fitted %.3f mm (offset %.2f)\n",
        mean(qfix), mean(qfree), mean(ctr_nom), mean(ctr_fit), mean(ctr_nom) - mean(ctr_fit))
@printf("  free/fixed ratio = %.4f +- %.4f  (constant => a pure SCALE error, not a shape error)\n",
        mean(Mfree ./ Mfix), std(Mfree ./ Mfix))

near(H, M, t) = M[argmin(abs.(H .- t))]
@printf("\nat 6.975 T:  fixed %.4f | free %.4f | 2.5 K VSM %.4f | digitised %.4f  (PRL SI ~1.65)\n",
        near(H04, Mfix, 6.975), near(H04, Mfree, 6.975), near(H25, M25, 6.975),
        near(Hdig, Mdig, 6.975))

fig = Figure(size = (1500, 1000))
Label(fig[0, 1:2],
      "MPMS3 centring correction: the 0.42 K M(H) scale problem was a choice of moment column";
      fontsize = 17, font = :bold)

# ---- 1  the curves ----------------------------------------------------------------
ax1 = Axis(fig[1, 1], xlabel = "field (T)", ylabel = "M (uB / Yb)",
           title = "Free centring agrees with an independent instrument; fixed does not")
scatter!(ax1, Hdig, Mdig; color = (:black, 0.55), markersize = 7,
         label = "digitised curve (used until now)")
lines!(ax1, H04, Mfix; color = :indianred, linewidth = 2.8,
       label = "0.42 K, FIXED centring")
lines!(ax1, H04, Mfree; color = :dodgerblue, linewidth = 3.0,
       label = "0.42 K, FREE centring (correct)")
lines!(ax1, H25, M25; color = :seagreen, linewidth = 2.4, linestyle = :dash,
       label = "2.5 K, DynaCool VSM (different crystal)")
hlines!(ax1, [3.436 / 2]; color = :grey40, linestyle = :dot, linewidth = 2)
text!(ax1, 0.2, 3.436/2 + 0.03; text = "g_par*S = 1.718 at g_par = 3.436 (PRL SI)",
      fontsize = 10, color = :grey30)
xlims!(ax1, 0, 7.3); axislegend(ax1; position = :rb, labelsize = 10)

# ---- 2  the ratio is constant, so it is a scale error -----------------------------
ax2 = Axis(fig[1, 2], xlabel = "field (T)", ylabel = "free / fixed",
           title = "The error is a CONSTANT factor: a scale problem, not a shape problem")
ok = (Mfix .> 0.02) .& isfinite.(Mfree)
r = Mfree[ok] ./ Mfix[ok]
scatter!(ax2, H04[ok], r; color = :purple, markersize = 7)
hlines!(ax2, [mean(r)]; color = :black, linestyle = :dash, linewidth = 2)
text!(ax2, 1.0, mean(r) + 0.004;
      text = @sprintf("mean = %.4f, spread %.4f", mean(r), std(r)), fontsize = 11)
ylims!(ax2, mean(r) - 0.06, mean(r) + 0.06); xlims!(ax2, 0, 7.3)
text!(ax2, 0.3, mean(r) - 0.05;
      text = "Because it is constant, the fitted amplitude A_M absorbed it entirely --\n" *
             "so no SHAPE fit in this repo was ever affected, but the absolute scale was.",
      fontsize = 10, color = :grey25)

# ---- 3  the file's own verdict ----------------------------------------------------
ax3 = Axis(fig[2, 1], xlabel = "field (T)", ylabel = "fit quality",
           title = "The file says which column to trust")
lines!(ax3, H04, qfix; color = :indianred, linewidth = 2.6, label = "DC Fixed Fit (poor)")
lines!(ax3, H04, qfree; color = :dodgerblue, linewidth = 2.6, label = "DC Free Fit (good)")
ylims!(ax3, 0, 1.08); xlims!(ax3, 0, 7.3)
axislegend(ax3; position = :rc, labelsize = 10)
ax3b = Axis(fig[2, 1], yaxisposition = :right, ylabel = "sample centre (mm)",
            ygridvisible = false)
hidespines!(ax3b); hidexdecorations!(ax3b); linkxaxes!(ax3, ax3b)
lines!(ax3b, H04, ctr_nom; color = :grey30, linewidth = 2, linestyle = :dot)
lines!(ax3b, H04, ctr_fit; color = :darkorange, linewidth = 2.4)
text!(ax3b, 1.6, mean(ctr_fit) + 0.9;
      text = @sprintf("nominal %.2f mm (dotted)\nfitted %.2f mm (orange)\noffset %.2f mm",
                      mean(ctr_nom), mean(ctr_fit), mean(ctr_nom) - mean(ctr_fit)),
      fontsize = 10, color = :darkorange)

# ---- 4  the physics sanity check --------------------------------------------------
ax4 = Axis(fig[2, 2], xlabel = "field (T)", ylabel = "M (uB / Yb)",
           title = "Colder must give MORE magnetisation -- only the corrected curve does")
scatter!(ax4, Hdig, Mdig; color = (:black, 0.5), markersize = 7, label = "digitised 0.42 K")
lines!(ax4, H04, Mfree; color = :dodgerblue, linewidth = 3.0, label = "0.42 K corrected")
lines!(ax4, H25, M25; color = :seagreen, linewidth = 2.4, linestyle = :dash, label = "2.5 K")
xlims!(ax4, 0, 3.2); ylims!(ax4, 0, 1.2)
axislegend(ax4; position = :rb, labelsize = 10)
text!(ax4, 0.15, 1.02;
      text = "Corrected 0.42 K sits ABOVE 2.5 K at low field and converges by ~7 T,\n" *
             "where Zeeman energy far exceeds kT. The digitised curve sat BELOW\n" *
             "2.5 K everywhere, which no paramagnet can do.",
      fontsize = 10, color = :grey25)

Label(fig[3, 1:2],
      "Also in this file: the `Temperature (K)` column reads 1.56 K -- that is the MPMS3 CHAMBER. " *
      "The sample temperature lives in a column named `He3 temp` and reads 0.42-0.50 K, matching " *
      "the filename. A reduction that used `Temperature (K)` would be wrong by a factor of four.";
      fontsize = 11, color = :grey30)

out = joinpath(FDIR, "mpms3_centering_correction.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
