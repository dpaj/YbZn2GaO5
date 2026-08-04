# PPMS magnetization, 2.5 K, 0-14 T

**Drop the two raw `.DAT` files in this directory without renaming them.** Original instrument
filenames are provenance; a `README` is the right place for description, not the filename.

Measured by **Aya Rutherford** on a **PPMS**, on **our** YbZn2GaO5 crystals — a *different sample*
from the one in the Bag et al. PRL. Two orientations, referred to here as "para" and "perp".

## Why this dataset exists (three purposes, in order of value)

1. **Discriminate g-anisotropy from exchange-anisotropy.** This is the strongest use and it was not
   the original motivation. Our model puts *all* anisotropy in the g tensor and is isotropic
   Heisenberg in exchange; the published model uses XXZ with `Delta ~ 1.35`, i.e. anisotropy in the
   *exchange*. Both give an anisotropic magnetic response, but with different field dependence, so
   the **ratio of the two orientations** is a direct test of which picture is right. Note the model
   cannot currently predict H perpendicular to c at all: `gperp_ratio` in
   `configs/best_fit_parameters.toml` is a *neutron matrix-element ratio*, not a magnetization
   g_perp, and `sunny_transverse_gxy = 1.0` is an intensity gauge. A real g_perp would have to be
   added.
2. **Sample-quality cross-check.** Reproduce the PRL Supplementary M(H) on a different crystal, to
   show our sample exhibits the same observables as the one behind the erroneous claim. This needs
   **no model at all** -- just an overplot -- so it is the cheapest result in the repo and should be
   done first. It requires the PRL SI curve to be digitized; that is not yet in this repo.
3. **Model validation, NOT a fit input.** Finite-temperature magnetization is expensive to compute
   and weakly sensitive to the parameters, so the target is semi-quantitative agreement. Do not
   optimise against it.

## What has to be recorded here before the data can be used quantitatively

- [ ] **Sample mass in mg, for each crystal separately.** Required, and nothing can be converted
      without it. PPMS reports moment in emu; the conversion is
      `mu_B per Yb = emu / (5585 * moles)` with `moles = mass_g / 453.53`, using
      `M(YbZn2GaO5) = 453.53 g/mol` (Yb 173.05 + 2xZn 130.76 + Ga 69.72 + 5xO 80.00) and one Yb per
      formula unit. `5585 emu/mol` is `N_A * mu_B`.
- [ ] **Which file is which orientation**, and what "para"/"perp" mean relative to the crystal
      axes. The model's easy axis is **c**, so presumably para = H parallel to c. Confirm rather
      than assume: it determines which file the existing model can be compared against at all.
- [ ] **Whether any normalisation or background subtraction is already applied** in the `.DAT`, or
      whether the moment column is raw instrument output. Sample-holder and diamagnetic
      contributions matter at 14 T.
- [ ] **Demagnetisation correction status**, and crystal shape if not corrected.

## A fourth benefit worth flagging

The existing low-temperature curve
(`data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv`) is **digitised from a figure** and
saturates near 1.12 uB/Yb, well below `gzz * S = 1.75` at `gzz = 3.5`. Whether that gap is a
normalisation error in the digitised data has been an open question in this repo from early on (see
CLAUDE.md open threads, the excess linear M(H) term).

**PPMS data with a known sample mass gives an absolute `mu_B` per Yb with no digitisation step**, so
it can settle that question independently. Different temperature, so not a direct substitute -- but
it does bound the absolute scale.

## Computational note

2.5 K is **not** the regime where temperature can be ignored. CLAUDE.md records that temperature is
irrelevant to M(H) *at 0.42 K* (T = 0 gives rms 0.0182 against 0.0176 classical), and that finding
must not be carried here: at 2.5 K, `kT = 0.215 meV` against an exchange scale of ~0.5-1 meV and a
Zeeman scale of 1-3 meV. So this needs the **finite-temperature Langevin path**
(`scripts/sunny_largecell_mvh_classical.jl`, with `relax_before_thermalize = true` -- omitting it
once gave M(7 T) ~ 0.02 instead of ~1.73), not the T = 0 minimiser that `sv_mvh_curve` uses.

Classical statistics is a *better* approximation at 2.5 K than at 0.07 K, so this is a more
favourable regime for the method than the low-temperature data, not a harder one.
