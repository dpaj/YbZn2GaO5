#!/usr/bin/env julia
# DOES TEMPERATURE PRODUCE THE MYSTERY UPWARD SLOPE IN M(H)?  ~30 s, no Sunny, no KPM.
#
#   julia --project=. scripts/check_thermal_slope_disorder.jl
#
# THE HYPOTHESIS, and it is a good one because it is a PREDICTION of the disorder model rather than
# a patch on it. If quenched disorder drives some local excitations close to hw = 0, then at finite
# temperature those are thermally populated and the moment is reduced. Raising the field gaps them
# out, so the reduction SHRINKS with field -- which adds a POSITIVE contribution to dM/dB. That is
# exactly the sign of the excess linear term the M(H) fit keeps wanting, and the term that the
# crystal field says is not Van Vleck (fitted 0.19 uB/T against 0.0171 from the CEF, ~11x).
#
# The mechanism is disorder-specific: with sigma_gzz = 0.80 on a mean of 3.5, the LOW-g tail has a
# much smaller Zeeman splitting than the mean site, so those sites depolarize first. A clean model
# has no such tail.
#
# WHY NOT LSWT + BOSE, which the repo already has. _lswt_magnon_depletion is only valid ABOVE
# saturation -- it needs a gapped field-polarized state, and below saturation sum n_B diverges. Strong
# gzz disorder pushes saturation well past the clean B_sat, so the ENTIRE low-field region, which is
# where this effect lives and where the model's absolute error is worst (ratio 1.78 at 0.25 T against
# 1.10 at 7 T), is out of its reach. Using it here would answer only where the answer is boring.
#
# SO: single-site QUANTUM statistics with mean-field exchange -- the minimal tool that is correct in
# the two limits that matter. It keeps the exact two-level Brillouin form for the single-site physics,
# which is what fails classically, and approximates only the exchange:
#
#     h_i    = g_i*muB*B - z*J1*<s>            (AFM exchange opposes the field)
#     s_i    = (1/2) tanh(h_i / (2 kT))        (exponential saturation, not algebraic)
#     <s>    = mean_i s_i                     (solved self-consistently)
#     M      = mean_i g_i * s_i
#
# VALIDATION IS BUILT IN, and this is not optional. A roll-our-own must contain a limit that reduces
# to something independently known, or a units slip goes undetected -- the DGX lost a factor of g^2
# exactly this way and caught it only because such a check existed. Here: at J1 = 0 and sigma_gzz = 0
# this must reduce to the S = 1/2 Brillouin function M = g*(1/2)*tanh(g*muB*B/(2kT)) to machine
# precision, and at T -> 0 it must give M = <g>/2 exactly.
#
# LIMITS, stated up front. Mean-field on a FRUSTRATED triangular lattice overestimates order and
# misses fluctuations, so this is quantitative where Zeeman beats exchange and only indicative below
# that -- the crossover is near 2 T. It also has no Van Vleck (that is a separate additive term) and
# no spatial correlation of the disorder. The question here is only whether the thermal mechanism has
# the right SIGN and MAGNITUDE to matter, which mean-field can answer.
#
# TEMPERATURES. The AC susceptibility is at 20 mK and 450 mK -- NOT 0.07 K, which is the neutron
# temperature; all three are included. Note 450 mK is essentially the DC 420 mK, so the AC pair
# BRACKETS the DC point and the 20 mK vs 450 mK comparison is the cleanest available temperature
# lever, being one instrument rather than two.

using Printf, Statistics, Random, CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const MUB = 0.05788381806     # meV/T
const KB  = 0.08617333262     # meV/K
const ZNN = 6                 # triangular lattice nearest neighbours
const FDIR = SV.sv_repo_path(REPO_ROOT, "results/figures/sunny_validation/thermal_slope")
const TDIR = SV.sv_repo_path(REPO_ROOT, "results/feature_tables/sunny_validation/thermal_slope")
mkpath(FDIR); mkpath(TDIR)

"Per-site g values: Gaussian, clipped at zero, matching sv_apply_disorder! conventions."
function g_sites(gbar, sigma_gzz; n=4_000, seed=20260805)
    sigma_gzz <= 0 && return fill(Float64(gbar), n)
    rng = MersenneTwister(seed)
    return [max(0.0, gbar + sigma_gzz * randn(rng)) for _ in 1:n]
end

"""
Self-consistent mean-field magnetization for effective S = 1/2 with g disorder.

`T_K = 0` is handled exactly rather than as a limit: every site with h_i > 0 is fully polarized,
so M = mean(g_i)/2 and no iteration is needed.
"""
function mvh_meanfield(B, gs, J1, T_K; tol=1e-11, maxit=400)
    # THE T = 0 BRANCH USED TO RETURN mean(gs)/2 UNCONDITIONALLY, AND THAT WAS A BUG that corrupted
    # the headline result of this script. It assumed every site fully polarized ALONG the field, but
    # the finite-T formula at T -> 0 gives s_i = -1/2 for any site with h_i < 0, i.e. g_i*muB*B <
    # z*J1*<s>. At 2.5 T with gbar = 3.5, sigma_gzz = 0.8, J1 = 0.15 that is roughly 19% of sites.
    # So the reference sat ABOVE the true T -> 0 limit of this very function, and the "thermal
    # deficit" measured against it absorbed that gap.
    #
    # The symptom was visible and I missed it: at 20 mK, kT = 0.0017 meV against a 0.5 meV Zeeman,
    # so thermal population is nil -- yet the "thermal" slope came out 0.122 uB/T. A thermal effect
    # cannot survive T -> 0. The DGX caught it by asking exactly that question of the 20 mK row.
    #
    # Now T = 0 is evaluated with the SAME self-consistency at a negligible temperature, so the
    # reference is the true limit of the model rather than an idealisation of it, and what remains in
    # the deficit is thermal by construction.
    T_eff = T_K <= 0 ? 1e-6 : T_K
    twokT = 2 * KB * T_eff
    s = 0.5
    it = 0
    for k in 1:maxit
        it = k
        snew = 0.0
        for g in gs
            h = g * MUB * B - ZNN * J1 * s
            snew += 0.5 * tanh(h / twokT)
        end
        snew /= length(gs)
        abs(snew - s) < tol && (s = snew; break)
        s = 0.5 * (s + snew)            # damped, so the AFM self-consistency cannot oscillate
    end
    M = mean(g * 0.5 * tanh((g * MUB * B - ZNN * J1 * s) / twokT) for g in gs)
    # VALIDITY. A single uniform <s> cannot represent an antiferromagnet: below the field where
    # Zeeman overcomes exchange there is no uniform polarized solution at all, the true state being
    # a canted/120-degree structure with only a small uniform component. Forced to iterate anyway,
    # the AFM term flips every spin and the self-consistency lands on an unphysical fixed point with
    # M < 0. So a point is trustworthy only where Zeeman genuinely beats exchange AND M > 0.
    zeeman   = mean(gs) * MUB * B
    exchange = ZNN * J1 * abs(s)
    return (; M, s, iters = it, valid = M > 0 && zeeman > exchange,
              zeeman, exchange)
end

# ---------------------------------------------------------------------------------
# VALIDATION FIRST. If these fail nothing below is worth reading.
# ---------------------------------------------------------------------------------
println("VALIDATION -- the two limits this must reproduce")
brillouin(g, B, T) = g * 0.5 * tanh(g * MUB * B / (2 * KB * T))
# Wrapped in a function: a bare top-level `for` writing to `maxerr` hits Julia soft scope and dies
# with UndefVarError. See the gotcha in CLAUDE.md.
function validate_brillouin()
    maxerr = 0.0
    for (g, B, T) in ((3.5, 1.0, 0.42), (3.0, 7.0, 2.5), (3.8, 0.3, 0.02), (2.0, 14.0, 1.0))
        mine = mvh_meanfield(B, fill(g, 1), 0.0, T).M
        ref  = brillouin(g, B, T)
        maxerr = max(maxerr, abs(mine - ref))
        @printf("  Zeeman only  g=%.1f B=%5.1f T=%.2f K :  %.12f vs Brillouin %.12f   diff %.2e\n",
                g, B, T, mine, ref, mine - ref)
    end
    return maxerr
end
maxerr = validate_brillouin()
@printf("  worst deviation from the Brillouin function: %.2e  ->  %s\n", maxerr,
        maxerr < 1e-12 ? "PASS" : "FAIL -- do not trust anything below")
# The T -> 0 check must be done at a field high enough that EVERY site has h_i > 0, or it is not a
# check of the limit but of the idealisation that used to be hard-coded here. At 20 T that holds; at
# 7 T it does not quite, and the residual is reported rather than asserted -- it is the size of the
# effect the old T = 0 branch was silently attributing to temperature.
gs_t = g_sites(3.5, 0.8)
m20 = mvh_meanfield(20.0, gs_t, 0.15, 0.0).M
@printf("  T -> 0 at 20 T (all h_i > 0): M = %.10f, <g>/2 = %.10f, diff %.2e  ->  %s\n",
        m20, mean(gs_t)/2, m20 - mean(gs_t)/2,
        abs(m20 - mean(gs_t)/2) < 1e-9 ? "PASS" : "FAIL")
for Bc in (2.5, 4.0, 7.0)
    r = mvh_meanfield(Bc, gs_t, 0.15, 0.0)
    frac = count(g -> g * MUB * Bc < ZNN * 0.15 * r.s, gs_t) / length(gs_t)
    @printf("  T -> 0 at %4.1f T: M = %.6f vs <g>/2 = %.6f  (%.1f%% of sites have h_i < 0)\n",
            Bc, r.M, mean(gs_t)/2, 100 * frac)
end
println()
maxerr < 1e-12 || error("Brillouin limit failed; the implementation is wrong.")

# ---------------------------------------------------------------------------------
const GBAR, SIGG, J1 = 3.50, 0.80, 0.15
const TEMPS = [(0.0, "T = 0"), (0.020, "20 mK (AC cold)"), (0.070, "70 mK (neutron)"),
               (0.420, "420 mK (DC MPMS3)"), (0.450, "450 mK (AC warm)"), (2.50, "2.5 K (DC VSM)")]
# STARTS AT 2.5 T, NOT AT ZERO, AND THAT LIMITATION IS THE MAIN RESULT OF THIS SCRIPT.
#
# A single uniform <s> cannot represent an antiferromagnet. Below the field where Zeeman overcomes
# exchange -- z*J1*<s> = 6*0.15*0.5 = 0.45 meV against gbar*muB*B = 0.04 meV at 0.2 T -- there is no
# uniform polarized solution at all; the true state is canted with only a small uniform component.
# Iterated anyway, the AFM term flips every spin and the self-consistency lands on an unphysical
# fixed point: an earlier version of this script ran from 0.05 T and reported a "thermal deficit" of
# 3.48 uB against a total moment of 1.74, i.e. M < 0. Those numbers were meaningless.
#
# NOTE WHAT THE BUILT-IN VALIDATION DID AND DID NOT CATCH. The Brillouin check passes to machine
# precision -- and could not possibly have caught this, because it sets J1 = 0, which is exactly the
# limit where the broken term vanishes. THE LESSON, and it generalises: a limit check must exercise
# the term you are worried about, or it certifies only the part you were not worried about.
const BS = collect(range(2.5, 7.0; length=100))
gs  = g_sites(GBAR, SIGG)
gcl = g_sites(GBAR, 0.0)

@printf("Model: gbar = %.2f, sigma_gzz = %.2f, J1 = %.3f meV, z = %d, %d sites\n",
        GBAR, SIGG, J1, ZNN, length(gs))
@printf("kT at 20 mK = %.5f meV, at 420 mK = %.5f meV;  gbar*muB*B at 0.2 T = %.5f, at 1 T = %.5f\n\n",
        KB*0.020, KB*0.420, GBAR*MUB*0.2, GBAR*MUB*1.0)

curves = Dict{String,Vector{Float64}}()
clean  = Dict{String,Vector{Float64}}()
for (T, lbl) in TEMPS
    curves[lbl] = [mvh_meanfield(B, gs,  J1, T).M for B in BS]
    clean[lbl]  = [mvh_meanfield(B, gcl, J1, T).M for B in BS]
end

# The thermal DEFICIT relative to T = 0, and the extra slope it induces.
println("THERMAL DEFICIT relative to T = 0, and the slope it adds (disordered model)")
@printf("  %6s", "B(T)")
for (_, lbl) in TEMPS[2:end]; @printf(" %14s", split(lbl)[1]); end
println()
for Bt in (2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0)
    j = argmin(abs.(BS .- Bt))
    @printf("  %6.2f", BS[j])
    for (_, lbl) in TEMPS[2:end]
        @printf(" %14.5f", curves[TEMPS[1][2]][j] - curves[lbl][j])
    end
    println()
end

"Extra dM/dB from the thermal deficit shrinking with field: -d(deficit)/dB."
function thermal_slope(lbl, blo, bhi)
    i = argmin(abs.(BS .- blo)); j = argmin(abs.(BS .- bhi))
    d(k) = curves[TEMPS[1][2]][k] - curves[lbl][k]
    return -(d(j) - d(i)) / (BS[j] - BS[i])
end
println("\nEXTRA SLOPE from the thermal mechanism, uB/T, over two windows")
println("  compare against the FITTED excess 0.186 uB/T and the crystal-field Van Vleck 0.0171")
@printf("  %20s %16s %16s\n", "temperature", "2.5-4 T", "4-7 T")
for (_, lbl) in TEMPS[2:end]
    @printf("  %20s %16.5f %16.5f\n", lbl, thermal_slope(lbl, 2.5, 4.0), thermal_slope(lbl, 4.0, 7.0))
end

println("\nIS THE EFFECT DISORDER-SPECIFIC? deficit at 420 mK, disordered vs clean")
@printf("  %6s %16s %16s %10s\n", "B(T)", "sigma_gzz=0.80", "sigma_gzz=0", "ratio")
for Bt in (2.5, 3.0, 4.0, 5.0, 6.0, 7.0)
    j = argmin(abs.(BS .- Bt))
    dd = curves[TEMPS[1][2]][j] - curves["420 mK (DC MPMS3)"][j]
    dc = clean[TEMPS[1][2]][j]  - clean["420 mK (DC MPMS3)"][j]
    @printf("  %6.2f %16.5f %16.5f %10s\n", BS[j], dd, dc,
            abs(dc) > 1e-9 ? @sprintf("%.2f", dd/dc) : "-")
end

println("\nTHE AC LEVER: 20 mK against 450 mK, one instrument, bracketing the DC 420 mK")
@printf("  %6s %14s %14s %14s\n", "B(T)", "M(20mK)", "M(450mK)", "difference")
for Bt in (2.5, 3.0, 4.0, 5.0, 6.0, 7.0)
    j = argmin(abs.(BS .- Bt))
    a = curves["20 mK (AC cold)"][j]; b = curves["450 mK (AC warm)"][j]
    @printf("  %6.2f %14.5f %14.5f %14.5f\n", BS[j], a, b, a - b)
end

open(joinpath(TDIR, "thermal_slope.csv"), "w") do io
    println(io, "B_T," * join(("M_" * replace(l, " " => "_") for (_, l) in TEMPS), ","))
    for j in eachindex(BS)
        @printf(io, "%.5f", BS[j])
        for (_, lbl) in TEMPS; @printf(io, ",%.8f", curves[lbl][j]); end
        println(io)
    end
end

# ---------------------------------------------------------------------------------
fig = Figure(size = (1500, 1050))
Label(fig[0, 1:2],
      "Can thermal population of disorder-induced low-lying levels explain the excess dM/dH?  " *
      "Mean-field, effective S = 1/2, quantum statistics.";
      fontsize = 15, font = :bold)

exp042 = SV.sv_read_magnetization_csv(
    SV.sv_repo_path(REPO_ROOT, "data/magnetization/YZGO_MvH_0p42K_Bparc_mpms3.csv"))

ax1 = Axis(fig[1, 1], xlabel = "field (T)", ylabel = "M (uB / Yb)",
           title = "M(H) versus temperature, sigma_gzz = 0.80")
cols = [:black, :navy, :teal, :crimson, :darkorange, :grey50]
for (k, (_, lbl)) in enumerate(TEMPS)
    lines!(ax1, BS, curves[lbl]; color = cols[k], linewidth = 2.3, label = lbl)
end
scatter!(ax1, exp042.B_T, exp042.M_muB_per_Yb; color = (:black, 0.35), markersize = 4,
         label = "measured 0.42 K")
axislegend(ax1; position = :rb, labelsize = 9); xlims!(ax1, 0, 7.1)

ax2 = Axis(fig[1, 2], xlabel = "field (T)", ylabel = "M(T=0) - M(T)  (uB / Yb)",
           title = "Thermal deficit: it SHRINKS with field, so it adds positive slope")
for (k, (_, lbl)) in enumerate(TEMPS[2:end])
    lines!(ax2, BS, curves[TEMPS[1][2]] .- curves[lbl]; color = cols[k+1], linewidth = 2.3,
           label = lbl)
end
hlines!(ax2, [0.0]; color = (:black, 0.4), linestyle = :dash)
axislegend(ax2; position = :rt, labelsize = 9); xlims!(ax2, 0, 7.1)

ax3 = Axis(fig[2, 1], xlabel = "field (T)", ylabel = "dM/dB  (uB / Yb / T)",
           title = "dM/dB -- is there a positive thermal contribution where it is needed?")
for (k, (_, lbl)) in enumerate(TEMPS)
    d = diff(curves[lbl]) ./ diff(BS)
    lines!(ax3, 0.5 .* (BS[1:end-1] .+ BS[2:end]), d; color = cols[k], linewidth = 2.2, label = lbl)
end
hlines!(ax3, [0.0171]; color = :green, linewidth = 2.0, linestyle = :dash,
        label = "CEF Van Vleck 0.0171")
hlines!(ax3, [0.186]; color = :red, linewidth = 2.0, linestyle = :dot,
        label = "fitted excess 0.186")
axislegend(ax3; position = :rt, labelsize = 8); xlims!(ax3, 0, 7.1)

ax4 = Axis(fig[2, 2], xlabel = "field (T)", ylabel = "M(T=0) - M(T)  (uB / Yb)",
           title = "Disorder-specific? 420 mK deficit, disordered vs clean")
lines!(ax4, BS, curves[TEMPS[1][2]] .- curves["420 mK (DC MPMS3)"]; color = :crimson,
       linewidth = 2.6, label = "sigma_gzz = 0.80")
lines!(ax4, BS, clean[TEMPS[1][2]] .- clean["420 mK (DC MPMS3)"]; color = :grey40,
       linewidth = 2.6, linestyle = :dash, label = "sigma_gzz = 0 (clean)")
axislegend(ax4; position = :rt, labelsize = 9); xlims!(ax4, 0, 7.1)

Label(fig[3, 1:2],
      "Mean-field on a FRUSTRATED lattice overestimates order and misses fluctuations, so this is " *
      "quantitative where Zeeman beats exchange (above ~2 T) and indicative below. Validated against " *
      "the S = 1/2 Brillouin function to machine precision at J1 = 0, and against M = <g>/2 at T = 0. " *
      "No Van Vleck here -- that is a separate additive term. The AC data are at 20 mK and 450 mK, " *
      "bracketing the DC 420 mK, which makes the 20-vs-450 mK pair the cleanest temperature lever " *
      "available since it is one instrument rather than two.";
      fontsize = 10, color = :grey30, word_wrap = true, tellwidth = false)

out = joinpath(FDIR, "thermal_slope_disorder.png")
save(out, fig; px_per_unit = 2)
println("\nwrote $out")
