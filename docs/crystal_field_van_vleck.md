# Crystal field and the Van Vleck term

```text
scripts/crystal_field_van_vleck.jl
julia --project=. scripts/crystal_field_van_vleck.jl
```

Purpose: decide from the published crystal field whether the linear high-field
term the M(H) fit wants can be single-ion Van Vleck at all, instead of leaving it
as a free phenomenological slope. Sunny's `stevens_matrices(7/2)` and
`spin_matrices(7/2)` supply the operators, so no operator algebra is hand-coded.

## Source

L. Zhao, T. Chen, M. B. Stone, Q. Zhang, C. L. Sarkis, S. M. Koohpayeh and
C. L. Broholm, *Quenched disorder in the triangular lattice antiferromagnet
YbZn2GaO5*, **Phys. Rev. B 113, 014437 (2026)**, DOI 10.1103/xn2m-1jb5, Table II.

This is YbZn2GaO5 itself, not an analogue. Stevens convention,
`H_CEF = sum B_n^m O_n^m`, D3d site symmetry allowing
(n,m) = (2,0), (4,0), (4,3), (6,0), (6,3), (6,6):

| B_n^m | value (meV) |
|---|---|
| B_2^0 | -0.78(3) |
| B_4^0 | 1.40(2)e-2 |
| B_4^3 | -8.2(2)e-1 |
| B_6^0 | 6.6(3)e-4 |
| B_6^3 | -3.0(5)e-2 |
| B_6^6 | 1.62(6)e-2 |

**Why this set is unusually trustworthy for a Van Vleck calculation.** Crystal
fields fit to level *energies* alone do not determine the eigenvectors, and Van
Vleck depends entirely on eigenvectors. This fit used INS peak positions **and
intensities**, plus the anisotropic saturation magnetization of Bag et al. as a
constraint — 7 observables for 6 parameters. The paper states the parameters
"cannot be uniquely determined without the inclusion of the Lande g-factor
constraints." That is exactly the condition that makes a CF-derived chi_VV usable.

## Validation

Reproducing the paper's own outputs from its parameters, as a check on both the
transcription and the Stevens convention:

| | calculated | published/measured |
|---|---|---|
| levels (meV) | 0, 38.5, 60.6, 95.5 | 0, 38.3, 60.6, 95.4 |
| g_par | 3.448 | 3.44 |
| g_perp | 3.190 | 3.04 |

Levels agree to 0.2 meV and g_par to 0.2%, so Sunny's Stevens normalization
matches the paper's. g_perp is 5% high, which is the one soft spot.

The paper's **point-charge** parameters, run through the same code, give levels
0, 38.2, 40.7, 91.1 meV and g_par = 1.14 — badly wrong. Useful as a sensitivity
check: the calculation is not insensitive to the parameters, so the agreement above
is meaningful.

## Result

    chi_VV^zz   = 0.01705 uB/T      (per Yb, H parallel c)
    chi_VV^perp = 0.00488 uB/T

Breakdown of the zz component by level: 38 meV contributes 0.01374, 95 meV
contributes 0.00331, and the 60 meV level contributes **exactly zero** because it
is the pure |±3/2> doublet and has no J_z matrix element with the ground doublet.

**Uncertainty.** Sampling the quoted B_n^m errors and keeping only draws that still
reproduce the measured levels within 4 meV and g_par within 5% (1237 of 4000 draws
accepted): `chi_VV^zz = 0.01705 +- 0.00067 uB/T`, range 0.0155 to 0.0189. So the
crystal-field determination is tight, about +-4%.

## The fitted linear term is not Van Vleck

The M(H) fit (`scripts/fit_mvh_only.jl`) wants a physical slope of
`A_M * chi_vv = 0.3712 * 0.0990 = 0.03675 uB/T`. Against the crystal field:

| | uB/T |
|---|---|
| chi_VV^zz from the crystal field | 0.0171 +- 0.0007 |
| linear term the M(H) fit wants | 0.0368 |
| **ratio** | **2.16** (1.94-2.37 over the accepted CF draws) |

**The fitted slope is about 2.2x larger than single-ion Van Vleck can supply**, and
the margin is far outside the crystal-field uncertainty. So the linear high-field
rise should not all be attributed to Van Vleck. Candidates for the excess:

- a normalization or background problem in the digitized M(H) curve (independently
  suspected on other grounds),
- a paramagnetic impurity contribution,
- model error: the disorder treatment may not capture the true approach to
  saturation, leaving a linear term to absorb the mismatch.

A sum-rule argument makes this harder to escape, not easier. Van Vleck is bounded by
`sum_excited |<n|Jz|0>|^2 = <Jz^2> - m_eff^2` with `m_eff = g_par/(2 g_J)`, so a
*larger* g_par leaves *less* room for Van Vleck. That matters because of the next point.

## An open tension worth chasing: is g_par really 3.44?

`g_par = 3.44` comes from saturation magnetization measured **to 14 T** (Bag et al.,
PRL 133, 266703 (2024)), and it enters the crystal-field fit as a constraint. But
Wu et al. (PRL 135, 046704 (2025)) measure to **45 T** in pulsed field and state
that M "begins to saturate at around 15 T", extrapolating back to a saturation
moment of **2.1(1) uB** after removing the Van Vleck contribution. That implies
`g_par = 2 * 2.1 = 4.2`, not 3.44.

Three things line up around 4.2 rather than 3.44:

- the pulsed-field saturation moment, 2.1(1) uB implies g_par ~ 4.2;
- the M(H) shape fit here prefers gzz = 4.67 (though it really constrains
  B_sat ~ 4.0 T, with gzz and J1 individually only half as well determined);
- the 14 T dataset used for g_par = 3.44 is, on Wu et al.'s own evidence, not
  saturated.

If g_par is nearer 4.2 the crystal-field fit would need redoing with the revised
constraint, its eigenvectors would change, and by the sum rule chi_VV would most
likely come out *smaller* — widening the discrepancy above rather than closing it.

## Disorder caveat

Zhao et al. refine Zn/Ga site mixing at x = 0.60(5), y = 0.35(5) and measure CF
levels broadened well beyond resolution (Lorentzian FWHM 5.9, 8.6, 8.4 meV against
2.4, 2.2, 2.0 resolution), which they attribute to a distribution of YbO6
geometries. They state explicitly that "site dependence of the Lande g factor for
Yb3+ is also anticipated."

So chi_VV is really a **distribution**, and the number above is its value for one
average environment. A 6 meV spread on a 38 meV level is ~15%, which propagates to
a comparable spread in chi_VV — still much smaller than the factor 2.2 discrepancy,
so it does not rescue the Van Vleck interpretation.

This also matters for the wider argument: the same paper that supplies these CF
parameters independently documents the site mixing and anticipates g-factor
randomness, which is the mechanism this project's disorder model assumes. By
contrast Bag et al. report "no observable intrinsic chemical site mixing" and
attribute the CF broadening to CEF-phonon (vibronic) coupling. That is a direct
literature disagreement on the central question.

## Related literature numbers

- YbZnGaO4 (sister compound, Ma et al., PRB 104, 224433 (2021)) needed Gaussian
  disorder widths of **Delta_J = 0.5 and Delta_g_par = 0.28** to reproduce spin-wave
  broadening — a quantitative precedent for the disorder widths used here
  (sigma_J = 0.5, sigma_gzz/gzz = 0.19).
- Bag et al. J1-J2 XXZ parameters: J2/J1 ~ 0.12, Delta ~ 1.35, J1 ~ 0.5 meV. Note
  this repo's model is isotropic Heisenberg with J1 ~ 0.24 meV, so both the
  anisotropy and the exchange scale differ from the published fit.
