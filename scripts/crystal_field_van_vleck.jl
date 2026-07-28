#!/usr/bin/env julia

# Crystal-field calculation for Yb3+ in YbZn2GaO5: level scheme, ground-doublet
# g-tensor, and Van Vleck susceptibility.
#
# Purpose: decide from the published crystal field whether the linear high-field
# term the M(H) fit wants (A_M * chi_vv = 0.0367 uB/T) can be single-ion Van Vleck
# at all, rather than leaving it as a free phenomenological slope.
#
# Source of parameters:
#   L. Zhao, T. Chen, M. B. Stone, Q. Zhang, C. L. Sarkis, S. M. Koohpayeh and
#   C. L. Broholm, "Quenched disorder in the triangular lattice antiferromagnet
#   YbZn2GaO5", Phys. Rev. B 113, 014437 (2026), DOI 10.1103/xn2m-1jb5, Table II.
#
# That fit is unusually well constrained for a crystal field: it used INS peak
# positions AND intensities, plus the anisotropic saturation magnetization of Bag
# et al. (PRL 133, 266703 (2024)) as a constraint, giving 7 observables for 6
# parameters. Crystal fields fit to energies alone do not determine the
# eigenvectors, and eigenvectors are exactly what Van Vleck depends on, so the
# g-factor constraint is what makes this set usable here.
#
# Stevens convention, H_CEF = sum_{n,m} B_n^m O_n^m, D3d site symmetry allowing
# n,m = (2,0), (4,0), (4,3), (6,0), (6,3), (6,6). Sunny supplies the Stevens
# operators, so no operator algebra is hand-coded here.
#
# Run with:
#   julia --project=. scripts/crystal_field_van_vleck.jl

using Printf
using LinearAlgebra
using Statistics
using Random
using Sunny

const J_YB = 7 / 2              # Yb3+ 2F7/2 ground multiplet
const G_J = 8 / 7               # Lande g factor
const MU_B = 0.05788381806      # meV / T
const KB = 0.08617333262        # meV / K

# --- published parameter sets ------------------------------------------------
# Keys are (n, m) with n the rank. Values in meV.
const CF_PUBLISHED = Dict(          # PRB 113, 014437 (2026), Table II, fitted
    (2, 0) => -0.78, (4, 0) => 1.40e-2, (4, 3) => -8.2e-1,
    (6, 0) => 6.6e-4, (6, 3) => -3.0e-2, (6, 6) => 1.62e-2)
const CF_PUBLISHED_ERR = Dict(      # quoted uncertainties, same table
    (2, 0) => 0.03, (4, 0) => 0.02e-2, (4, 3) => 0.2e-1,
    (6, 0) => 0.3e-4, (6, 3) => 0.5e-2, (6, 6) => 0.06e-2)
const CF_POINT_CHARGE = Dict(       # same table, point-charge model
    (2, 0) => 0.47, (4, 0) => 2.35e-2, (4, 3) => -9.2e-1,
    (6, 0) => 3.7e-4, (6, 3) => -2.13e-3, (6, 6) => 3.66e-3)
const CF_ARXIV = Dict(              # arXiv:2507.12592v1, pre-referee, DIFFERENT
    (2, 0) => -0.91, (4, 0) => 1.46e-2, (4, 3) => -7.5e-1,
    (6, 0) => 6.0e-4, (6, 3) => -3.1e-2, (6, 6) => 1.83e-2)

# Measured values to validate against.
const E_MEASURED = (38.3, 60.6, 95.4)   # meV, PRB Table IV
const G_PARA_LIT = 3.44                 # Bag et al., saturation magnetization to 14 T
const G_PERP_LIT = 3.04

# What the M(H) fit wants, as an actual dM/dB present in the data.
const SLOPE_FITTED = 0.3712 * 0.0990    # A_M * chi_vv, uB/T

# -----------------------------------------------------------------------------

"Crystal-field Hamiltonian from Stevens parameters, in the |J, m> basis."
function cf_hamiltonian(B::Dict; J=J_YB)
    O = stevens_matrices(J)
    N = Int(2J + 1)
    H = zeros(ComplexF64, N, N)
    for ((n, m), b) in B
        H .+= b .* O[n, m]
    end
    return Hermitian(H)
end

"""
Level scheme, ground-doublet g-tensor, and Van Vleck susceptibility.

Van Vleck is the second-order Zeeman response,

    chi_VV^aa = 2 g_J^2 mu_B  sum_{n not in ground doublet} |<n|J_a|0>|^2 / (E_n - E_0)

averaged over the two ground states. Returned in uB/T per Yb. Only states outside
the ground doublet contribute: within-doublet matrix elements are first-order and
belong to the g-tensor instead.
"""
function cf_analyze(B::Dict; J=J_YB, doublet_tol=1e-6)
    H = cf_hamiltonian(B; J)
    F = eigen(H)
    E = real.(F.values) .- minimum(real.(F.values))
    V = F.vectors
    Jx, Jy, Jz = spin_matrices(J)

    # Group into Kramers doublets.
    groups = Vector{Vector{Int}}()
    for i in eachindex(E)
        placed = false
        for g in groups
            if abs(E[i] - E[g[1]]) < max(doublet_tol, 1e-6 * max(1.0, maximum(E)))
                push!(g, i); placed = true; break
            end
        end
        placed || push!(groups, [i])
    end
    sort!(groups; by=g -> E[g[1]])
    levels = [mean(E[g]) for g in groups]
    ground = groups[1]
    length(ground) == 2 || @warn "Ground level is not a doublet" n=length(ground)

    # Rotate the ground doublet so Jz is diagonal within it, giving the standard
    # |+>, |-> basis: then g_par = 2 g_J <+|Jz|+> and g_perp = 2 g_J |<+|Jx|->|.
    Vg = V[:, ground]
    jz_blk = Hermitian(Vg' * Jz * Vg)
    fg = eigen(jz_blk)
    Vg = Vg * fg.vectors
    up = argmax(real.(fg.values))
    dn = 3 - up
    vp, vm = Vg[:, up], Vg[:, dn]

    m_eff = real(vp' * Jz * vp)
    g_par = 2 * G_J * abs(m_eff)
    g_perp = 2 * G_J * abs(vp' * Jx * vm)

    # Van Vleck: sum over states outside the ground doublet.
    excited = setdiff(eachindex(E), ground)
    chi = Dict{Symbol,Float64}()
    contrib = NamedTuple[]
    for (sym, Jop) in ((:zz, Jz), (:xx, Jx), (:yy, Jy))
        tot = 0.0
        for gs in (vp, vm)
            for n in excited
                Δ = E[n] - E[ground[1]]
                Δ > 1e-9 || continue
                tot += abs2(V[:, n]' * Jop * gs) / Δ
            end
        end
        chi[sym] = 2 * G_J^2 * MU_B * tot / 2      # /2 averages the two ground states
    end
    # Per-level breakdown for the zz component.
    for (li, g) in enumerate(groups[2:end])
        s = 0.0
        for gs in (vp, vm), n in g
            s += abs2(V[:, n]' * Jz * gs)
        end
        s /= 2
        push!(contrib, (; level=li, E=levels[li + 1], sum_Jz2=s,
                          chi_zz=2 * G_J^2 * MU_B * s / levels[li + 1]))
    end

    return (; levels, E, V, g_par, g_perp, m_eff,
              chi_zz=chi[:zz], chi_perp=(chi[:xx] + chi[:yy]) / 2, contrib,
              ground_weights=[abs2(vp[i]) for i in 1:Int(2J + 1)])
end

function report(name, B; validate=false)
    r = cf_analyze(B)
    @printf("\n--- %s\n", name)
    @printf("    levels (meV)     : %s\n",
            join([@sprintf("%.1f", e) for e in r.levels], ", "))
    if validate
        @printf("    measured (meV)   : 0.0, %.1f, %.1f, %.1f\n", E_MEASURED...)
        d = [r.levels[i + 1] - E_MEASURED[i] for i in 1:min(3, length(r.levels) - 1)]
        @printf("    difference       : %s\n", join([@sprintf("%+.1f", x) for x in d], ", "))
    end
    @printf("    g_par, g_perp    : %.3f, %.3f", r.g_par, r.g_perp)
    validate && @printf("   (literature %.2f, %.2f)", G_PARA_LIT, G_PERP_LIT)
    println()
    @printf("    chi_VV zz, perp  : %.5f, %.5f uB/T\n", r.chi_zz, r.chi_perp)
    @printf("    zz breakdown     : %s\n",
            join([@sprintf("%.0f meV -> %.5f", c.E, c.chi_zz) for c in r.contrib], ", "))
    return r
end

function main()
    println("Crystal field for Yb3+ in YbZn2GaO5")
    println("Stevens convention, J = 7/2, g_J = 8/7, D3d site symmetry")
    @printf("Fitted M(H) linear term to explain: A_M*chi_vv = %.5f uB/T\n", SLOPE_FITTED)

    r = report("PRB 113, 014437 (2026) fitted [PRIMARY]", CF_PUBLISHED; validate=true)
    report("arXiv:2507.12592v1 fitted (pre-referee)", CF_ARXIV; validate=true)
    report("point-charge model (same paper)", CF_POINT_CHARGE; validate=true)

    println("\n================ does Van Vleck explain the fitted slope? ================")
    ratio = SLOPE_FITTED / r.chi_zz
    @printf("  chi_VV^zz from the published crystal field : %.5f uB/T\n", r.chi_zz)
    @printf("  linear term the M(H) fit wants             : %.5f uB/T\n", SLOPE_FITTED)
    @printf("  ratio fitted / crystal-field              : %.2f\n", ratio)
    println()
    if ratio > 1.5
        @printf("  The fitted slope is %.1fx larger than single-ion Van Vleck can supply.\n", ratio)
        println("  Something else contributes the high-field linear rise: a paramagnetic")
        println("  impurity, an unsubtracted background, or a normalization error in the")
        println("  digitized data. It should NOT all be attributed to Van Vleck.")
    elseif ratio > 0.7
        println("  Consistent within the uncertainties: Van Vleck can plausibly account")
        println("  for the fitted slope, and chi_vv can be fixed rather than fitted.")
    else
        println("  The fitted slope is SMALLER than the crystal field predicts, which")
        println("  would mean Van Vleck is being partly absorbed elsewhere in the model.")
    end

    # --- how much of this survives the underdetermination of the CF fit? -----
    # Sample the quoted parameter uncertainties and keep only draws that still
    # reproduce the measured level energies and g_par, then quote the spread.
    println("\n================ uncertainty from the crystal-field fit ================")
    println("  Sampling the quoted B_n^m uncertainties, keeping draws that still match")
    println("  the measured levels (within their errors) and g_par within 5%.")
    rng = MersenneTwister(20260728)
    nsamp, kept = 4000, NamedTuple[]
    for _ in 1:nsamp
        Bs = Dict(k => CF_PUBLISHED[k] + CF_PUBLISHED_ERR[k] * randn(rng) for k in keys(CF_PUBLISHED))
        rr = try
            cf_analyze(Bs)
        catch
            continue
        end
        length(rr.levels) >= 4 || continue
        ok = all(abs(rr.levels[i + 1] - E_MEASURED[i]) < 4.0 for i in 1:3) &&
             abs(rr.g_par - G_PARA_LIT) / G_PARA_LIT < 0.05
        ok && push!(kept, (; chi_zz=rr.chi_zz, chi_perp=rr.chi_perp,
                             g_par=rr.g_par, E1=rr.levels[2]))
    end
    if length(kept) < 10
        @printf("  only %d of %d draws satisfied the constraints — widen them to say more\n",
                length(kept), nsamp)
    else
        cz = [k.chi_zz for k in kept]
        @printf("  %d of %d draws accepted\n", length(kept), nsamp)
        @printf("  chi_VV^zz = %.5f +- %.5f uB/T   (range %.5f - %.5f)\n",
                mean(cz), std(cz), minimum(cz), maximum(cz))
        @printf("  ratio fitted/CF over the accepted set: %.2f - %.2f\n",
                SLOPE_FITTED / maximum(cz), SLOPE_FITTED / minimum(cz))
        @printf("  g_par over accepted set: %.3f +- %.3f\n",
                mean(k.g_par for k in kept), std([k.g_par for k in kept]))
    end

    println("\n================ notes and caveats ================")
    println("  * Zhao et al. refine Zn/Ga site mixing at x = 0.60(5), y = 0.35(5) and")
    println("    measure CF levels broadened well beyond resolution (Lorentzian FWHM")
    println("    5.9, 8.6, 8.4 meV against 2.4, 2.2, 2.0 resolution). They explicitly")
    println("    anticipate a site-dependent Lande g factor. So chi_VV is really a")
    println("    DISTRIBUTION and this is its mean over one average environment.")
    println("  * A ~6 meV spread on a 38 meV level is ~15%, which propagates to a")
    println("    comparable spread in chi_VV, well below the discrepancy above.")
    println("  * g_par = 3.44 comes from saturation magnetization to 14 T. Wu et al.")
    println("    (PRL 135, 046704 (2025)) find with pulsed field to 45 T that M does")
    println("    not saturate until ~15 T, and extrapolate 2.1(1) uB after removing")
    println("    Van Vleck, which implies g_par ~ 4.2 rather than 3.44. That is close")
    println("    to the 4.67 the M(H) shape fit prefers, and it matters here because")
    println("    g_par enters the CF fit as a constraint.")
    return r
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
