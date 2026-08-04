# YbZn2GaO5 (YZGO) — orientation

Read this first. It is the handoff document: it records the current scientific
state, the conventions that are easy to get wrong, and the open threads. The
per-topic detail lives in `docs/`.

## The scientific claim

YbZn2GaO5 was published as a **Dirac quantum spin liquid** (Bag et al., PRL 133,
266703 (2024); Wu et al., PRL 135, 046704 (2025)). The argument of this repository
is that the data are instead explained by **quenched disorder** from Zn/Ga site
mixing. Zhao et al. (PRB 113, 014437 (2026)) independently support the disorder
picture, refining site mixing at x = 0.60(5), y = 0.35(5) and explicitly
anticipating a distribution of Lande g factors. The two camps disagree head-on
about whether the crystal-field broadening is site mixing or CEF-phonon coupling.

**Who the argument is with.** In this repository "published" means **Bag et al., the
Haravifard group**. Work from the **Broholm group** identifying possible disorder in this
system is **aligned with the argument here, not a target of it** -- do not lump the two
together when writing anything up. (I have not verified which of the citations above maps to
which group, so check before attributing.)

### What is already robust, independent of the open caveats

Worth keeping in view, because the day-to-day work is currently deep in background systematics
and parameter precision, and neither of those threatens the thesis:

1. **YZGO is not well described by the published Hamiltonian.** Established, not pending.
2. **The field-saturated phase cannot be described by a single resolution-limited mode.**
   Disorder broadening is *required*. The clean control at the same exchange and g magnitudes
   gives a resolution-limited magnon against measured peaks 5-10x broader
   (`chi2_red` 378.7 against 24.1 on the 1D cuts, a thin line against a broad band in 2D), and
   the measured widths survive a move to the better-resolution Ei = 3.32 setting -- so the
   broadening is inhomogeneous, not instrumental.
3. **M(H) and the other observables follow the disorder model rather than the published one.**
   No plateau where the published parameters demand one, high-field diffraction showing only
   k = 0, non-Lorentzian magnon widths, and crystal field plus diffraction both requiring
   disorder to fit at all.

**These conclusions survive a systematic background bias and survive imprecision in the
parameter values.** The background work and the parameter refinement improve the quantitative
statement; they are not load-bearing for the qualitative one. Do not let a difficult background
region or an unresolved gzz digit get reported as though the thesis were in doubt.

Status: a minimal single-disordered-phase model reproduces M(H) semi-quantitatively
(rms 0.004 uB against a 1.12 uB signal) and the 1D neutron cuts qualitatively.

### Why the published model is decisively wrong, not merely a worse fit

This is the core of the argument and is **not** a matter of comparing chi2 values. Five
independent observations, four of them qualitative:

1. **High-field neutron diffraction sees only k = 0.** Field-induced on-site moments and
   **no k != 0 peaks whatsoever.** The published clean model requires structured non-k=0
   order; those phases are simply absent. **This dataset is not yet in the repo** -- it
   lives only in this note and in the experimenter's records, so do not assume a file
   exists for it.
2. **M(H) shows no plateau.** Simulations at the published parameters give a rich
   structured phase diagram with a clear magnetization plateau. The measurement is a
   monotonic, slope-changing approach to saturation -- exactly what the disorder model
   produces and nothing like the clean model.

   Provenance worth knowing, because it is not obvious from the repo. **Our M(H) is the only
   available low-temperature set**: the original publication either did not measure
   low-temperature magnetization or did not publish it. The same group later presented
   magnetization at an APS Global Summit talk and explained it with a *phenomenological
   singlet model*, which the experimenter here judges to be unfounded -- and critically,
   **their M(H) curves are highly similar to ours**, so the measurement itself is
   reproducible across groups even though the interpretation is disputed.

   The published Hamiltonian's failure on M(H) is **far from even qualitative** agreement and
   has already been demonstrated with separate scripts held outside this repo. It does NOT
   need re-deriving here, and there is no need to build a published-Hamiltonian M(H)
   comparison in this codebase.
3. **The Bag et al. zero-field neutron fits are not even qualitatively right.**
4. **The high-field magnon widths exceed the instrumental resolution by a lot, and are
   NOT Lorentzian.** So the broadening is inhomogeneous -- a distribution of parameters --
   not lifetime/time broadening. This is direct evidence for disorder as the mechanism.
5. **Crystal field and diffraction data must use a disorder model to fit at all.**

Consequences for how to work here. XXZ anisotropy (the literature uses Delta ~ 1.35) may
eventually prove necessary, but it is **unlikely to be the main story, and if needed it
goes ON TOP OF the disorder model, not instead of it.** Do not read a poor chi2 on the
dispersive cuts as support for the published model -- see the K/M misfit note under
Established results.

Also note that **J1 and most other parameters were inherited from the earlier analytical
model**, which cannot represent mode mixing between the different exchange environments.
Sunny/KPM can. So the mean of the J1 distribution shifting once a proper disorder
distribution is used is *expected*, not a red flag. The old analytical co-refinement had a
known unresolved discrepancy of exactly this kind: the isolated "flat" neutron mode was
taking up work that mode mixing of the exchange-broadened spectrum should have done.

## Layout

```text
configs/   TOML run controls, one per script; deep-merged onto a base config
data/      experimental inputs (tracked; data/raw and data/generated are ignored)
docs/      the real documentation - see the reading order below
results/   generated figures and tables — GITIGNORED, all regenerable
scripts/   runnable entry points; scripts/dev holds benchmarks; scripts/legacy the old co-fit
src/       sunny_validation.jl (large), feature_extraction.jl, parameters.jl
test/      runtests.jl — physics-invariant regression tests (see Housekeeping)
```

`../references/` (a sibling of this repo, **not** tracked) holds the five key
papers. A clone will not have them.

Reading order for the current work: `docs/largecell_mvh_classical.md`,
`docs/crystal_field_van_vleck.md`, `docs/benchmarking.md`, then
`docs/companion/` for the older analytical layer.

## Working principle: never bury an assumption

Where a choice has to be made and cannot be measured -- a window edge, an interpolation form, a
weighting -- **make the choice explicit, run the defensible alternatives, and carry the spread
forward as an uncertainty.** It usually costs little and it makes the assumption auditable
instead of invisible.

The corollary matters as much: an uncertainty band built by varying parameters *within* an
assumption does NOT cover the assumption being wrong. Stage 1's 18-variant envelope spans window
edges and interpolation form; it says nothing about whether "no known process peaks at 1.5 meV"
is true. State that limit whenever the band is quoted.

The failure mode this replaces is real and recurred several times here: fixing a bias by editing
the objective, which trades a measured bias for an unmeasured one and hides the choice inside a
number.

## Conventions that are easy to get wrong

- **Effective crystal, not the CIF.** `sv_effective_triangle_crystal` builds a
  **P1** one-site triangular net. Because it is P1, Sunny will not generate the
  triangular shells from one representative bond — the J1 and J2 offsets are
  hard-coded in `sv_j1_shell_offsets` / `sv_j2_shell_offsets`. The test suite pins
  them against the analytical J(q) form factors.
- **Isotropic Heisenberg.** J1 and J2 are scalar exchanges. There is no Delta / XXZ
  anisotropy and no single-ion anisotropy anywhere; all anisotropy lives in the
  g-tensor. Note the published fits use an XXZ model with Delta ~ 1.35, so this
  model differs from the literature.
- **`sunny_transverse_gxy = 1.0` is an intensity gauge, not a physical gperp.**
  `ssf_perp` needs nonzero transverse moment components or the inelastic intensity
  vanishes numerically. Never report it as a refined g-factor.
- **Sign.** Sunny's moment is `-g S`, so `sv_m_parallel_uB_per_site` returns the
  calculator sign and `[largecell].moment_sign = -1.0` converts to the
  experimental convention. Always apply it.
- **Disorder.** `sigma_J` is **fractional** (`J*(1 + sigma_J*randn)`), `sigma_gzz`
  is **absolute** and clipped at zero. `sv_apply_disorder!` takes a `realization`
  keyword; `realization = 0` reproduces the original single-realization seed
  bit-for-bit, so published neutron/KPM results are unchanged.
- **Two magnetization normalizations exist.** The Sunny path divides by `(1+r2)`;
  the legacy analytical co-fit does not. `magnetization_global_scale = 0.4189` came
  from the un-normalized fit, so `A_M` is not comparable between them.
- **`[largecell]` in the base config is superseded** except for `moment_sign` and
  the `include_*_disorder` flags. Its `4x4x1` drives only the legacy script.

## The data, and how far to trust it

**Read this before interpreting any neutron `chi2` in this repo.** The background under the
1D cuts is largely *constructed* rather than measured, and it is constructed precisely where
the interesting signal lives.

`sv_min_over_fields_background_raw` takes the minimum over 0, 9 and 14 T in **two anchor
windows only** -- `[0, 0.75]` meV and `E > 2.5` meV -- and PCHIP-interpolates everything
between. The neutron fit window is `[0.5, 3.0]` meV, so **1.75 of its 2.5 meV, roughly 70% of
the fit window, rests on interpolated background.** There is also a known sharp magnet
background feature in the Ei = 4.65 meV data which required a bespoke removal scheme, and it
lives in that same gap.

The three cuts are also **not** corrected identically. Structured-residual windows exist for
`0p33_0p33_0` (1.675-2.375 meV) and `0p5_0_0` (1.825-2.425 meV), but there is **none for
`0_1_0`**, whose gap is therefore pure interpolation.

Consequences that should govern priorities here:

- **Sub-percent numerical convergence sits far below the systematic floor.** q-sampling error
  is under 1%, the realization floor is 12-15%, and the background across most of the fit
  window is an interpolation. Do not spend compute driving q error below a few tenths of a
  percent -- it cannot change a conclusion.
- **The residual that remains at 9 T -- data carrying weight around 1.5-2.4 meV that no
  parameter set reproduces -- sits dead centre in the unanchored gap**, and leftover magnet
  background would appear with exactly that sign. It may therefore be an artefact rather than
  missing physics. **Do not add model complexity (non-Gaussian disorder, XXZ) to chase it**
  without independent evidence that the feature is real.
- Higher-order experimental effects -- absorption, sample centring -- can shift the relative
  scale between one momentum or energy position and another. Cut-to-cut *amplitude* agreement
  is therefore weaker evidence than *lineshape* agreement.
- The arbiter is the **Ei = 3.32 meV** data, used as an approximate cross-check where
  kinematically accessible. Note that `(0,1,0)` is NOT accessible at that Ei, so the one cut
  with no structured-residual correction is also the one that cannot be cross-checked.

## Magnetization data: the observables, and one trap that cost a year

`data/magnetization/` and `data/ac_susceptibility/` now hold primary measurements. All
crystals come from the **same Haidong Zhou group growth** (measured by Aya Rutherford) as the
neutron crystal, so cross-comparison between them is legitimate.

| directory | instrument | T | field | orientation |
|---|---|---|---|---|
| `mpms3_0p4K/` | MPMS3 + He-3 insert | 0.42-0.50 K | 7 T | B \|\| c, 12.2 mg |
| `ppms_2p5K/` | DynaCool VSM | 2.5 K | 14 T | B \|\| c (4.81 mg), B ⟂ c (7.7 mg) |
| `ac_susceptibility/nhmfl/` | SCM1, dilution | ~29 mK | 18 T | both, 48 scans |

### THE MPMS3 CENTRING TRAP — read before reducing any MPMS3 file

**Every absolute magnetization number in this repo before 2026-08-04 was 1.49x too low.** The
cause is a data-reduction choice, not a measurement failure, and it does **not** affect any
shape fit — `A_M` is profiled out and the error is a constant — but it invalidates any absolute
value.

`Moment (emu)` is **EMPTY** in the MPMS3 file. The moment lives in two other columns, from two
different fits to the SQUID response:

- `DC Moment Fixed Ctr` — sample centre held at the nominal position
- `DC Moment Free Ctr` — centre floated as a fit parameter ← **USE THIS ONE**

The sample sat **3.06 mm** off centre (nominal 39.686 mm, fitted 36.623 mm), so the fixed fit
solved at the wrong position and under-read by **1.486 ± 0.004**. The file states its own
verdict: `DC Fixed Fit` quality **0.268** against `DC Free Fit` **0.957**.

Corrected, three independent numbers agree at ~7 T to four digits: MPMS3 0.42 K free centring
**1.6512**, DynaCool VSM 2.5 K on a *different crystal* **1.6522**, Bag et al. Supplement
**~1.65**.

**SECOND TRAP IN THE SAME FILE.** `Temperature (K)` reads **1.56 K** — that is the MPMS3
*chamber*. Sample temperature is in a column named **`He3 temp`**, reading 0.423-0.495 K. A
reduction using `Temperature (K)` is wrong by a factor of four. (The 0.42 K the M(H) protocol
has always assumed is correct.)

The physics check that settles it, and which was visible all along: the corrected 0.42 K curve
sits **above** 2.5 K at low field and converges by ~7 T where Zeeman ≫ kT. The old digitized
curve sat *below* 2.5 K at every field, which no paramagnet can do.

`data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv` is **superseded** — it was
digitized from a figure *and* came from the wrong column. It reproduces
`DC Moment Fixed Ctr` to 4 decimal places.

## Established results — do not re-derive

- **M(H) constrains B_sat, not J1 and gzz separately.** With `A_M` profiled out,
  gzz enters only through `B_sat = S*D_max/(gzz*mu_B)`; its amplitude role is
  absorbed. The ratio is ~2x better determined than either parameter. M(H) wants
  `B_sat ~ 4.0 T` where the by-eye neutron parameters give 5.1 T.
- **M(H) cannot constrain `sigma_J` at all.** Over the whole range 0 to 1 the rms
  moves by 2.7x the reproducibility floor. `sigma_J` must come from the spectra.
  `sigma_gzz` by contrast is tightly constrained. The two observables are therefore
  complementary, which is ideal for a co-optimization.
- **Temperature is irrelevant to M(H) at 0.42 K.** T = 0 and classical 0.42 K give
  rms 0.0182 vs 0.0176. Disorder, not temperature, produces both the rounding of
  the saturation and the extra initial slope, and both are present at T = 0.
- **12x12x1 is sufficient for M(H).** Realization scatter goes as N^-0.38 while
  cost goes as N^1.40, so accuracy is bought with realizations, not cell size
  (3.2x cheaper per unit accuracy). The ground-state texture has a correlation
  length of order one lattice constant, so the box is not limiting.
- **KPM is memory-bandwidth bound, but the ceiling is MACHINE-SPECIFIC.** On the
  Windows box CPU q-threading saturates near 3x; on the DGX it reaches 16.6x at 81
  chunks and is still rising, because a DGX has far more aggregate memory bandwidth.
  Both machines land at a similar *absolute* throughput (~0.05-0.06 s/q), which is
  the bandwidth wall; the speedup *factor* differs only because the serial baselines
  do. An earlier claim here that "the GPU port is the only effective lever" was
  wrong for the DGX: a single A100 in Float64 (0.065 s/q) actually LOSES to
  well-threaded DGX CPU (0.0495 s/q). **Multiple GPUs are the real lever** -- 4
  concurrent A100s give 3.91x aggregate, 97.8% of ideal. Host threading on top of
  one GPU helps ~1.2x at 36x36x1 but is net contention (0.88x) at 12x12x1, so it is
  system-size dependent.
- **Ground-state cost is small; an earlier claim that it was ~76% of runtime was a
  benchmarking error.** That figure was first-call JIT, not compute. Warm it is
  0.010 s at 12x12x1 and 0.18 s at 36x36x1 / 14 T (cold/warm ratio 579x at
  12x12x1), and warm times do scale with system size. There is no Amdahl cap on the
  GPU gain.
- **But `maxiters` is badly mis-set for 9 T.** At 9 T on 36x36x1 the energy is
  converged by ~1000 iterations, yet the minimizer never satisfies its convergence
  test (gradient plateaus near 8e-8 and is not even monotonic), so it burns the cap.
  E/site is identical to 8 decimals for maxiters of 1000, 5000, 20000 and 50000,
  while the time goes 0.89, 4.42, 17.8, 47.1 s. `maxiters = 50000` therefore wastes
  ~46 of 47 seconds at exactly the field where the neutron data lives. 14 T
  converges properly in 134 iterations; 9 T is pathological because the disordered
  system is only partially saturated there. **Judge convergence by E/site, not by
  the returned flag**, and do not raise `maxiters` to chase the flag.
- **`sigma_gzz`, not `sigma_J`, dominates both KPM cost and realization scatter.**
  Confirmed from two orthogonal scans on two machines: at fixed `sigma_J = 0.5`, cost rises
  **2.82x** across `sigma_gzz` 0 -> 1.0 (Windows); at fixed `sigma_gzz = 0.8`, it rises only
  **1.12x** across `sigma_J` 0 -> 1.0 (DGX). The realization floor behaves the same way --
  12.4% shape spread rms at `sigma_J = 0`, rising only ~19% by `sigma_J = 1.0`. The reason is
  a scale comparison: the Zeeman term is 1.7 meV at 9 T and 3.1 meV at 14 T and `sigma_gzz =
  0.8` spreads it by +-0.65 meV at 14 T, whereas the whole exchange bandwidth
  `S*[J(0)-J(q)]` is only ~0.5-1 meV at `J1 = 0.25`, `S = 1/2`. Both of us initially
  predicted 2-3x from `sigma_J`; the magnitude was right and the cause wrong.
  This is physics, not just performance: the g-factor distribution is the dominant disorder
  effect in the spectra, which is the mechanism Zhao et al. anticipated.
- **The neutron realization floor is 12-15%** (shape spread rms, per-cut range 0.080-0.191,
  6 realizations at 36x36x1). It is measured in SPECTRUM units, which is the wrong unit for a
  fit -- the number a fit consumes is the scatter in `chi2_red`, and that is still pending.
  Practical consequence already observed: at 4 realizations the Gamma surface separates
  `chi2_red` 9.4 from 50 trivially but does NOT resolve 9.42 from 9.60.
- **The Gamma cut alone prefers `gzz` = 3.30-3.45 with `sigma_gzz` = 0.80.** Complete 6x6
  surface on `(0,1,0)` at 9 and 14 T, where the exchange cancels identically so J does not
  enter to leading order. The by-eye `gzz = 3.8` scores 21.87 against 9.42, i.e. **2.3x
  worse**, and the published 3.44 sits in the minimum. `sigma_gzz = 0.8` is a genuine
  interior minimum, so the by-eye `sigma_gzz` was right and it is `gzz` that was wrong.
  Caveat: a mean `gzz` of 3.3 with `sigma_gzz = 0.8` is a 24% spread, so the convolved
  lineshape peak is NOT at `gzz*mu_B*B` -- do not compare the fitted mean directly against a
  peak-position estimate.
- **A staged factorization can be globally worse than its own starting point.** `gzz` shifts
  the mode energy at every q, so a `gzz` chosen from Gamma alone cannot see the cost it
  imposes at K and M; a 12x12x1 test landed 73% worse on the six-cut `chi2` while the
  "did Gamma move?" consistency check reported success. Always finish with a joint step, and
  note that an UNWEIGHTED six-cut sum is dominated by whichever cuts fit worst -- the DGX
  measured the same parameter point scoring ~9 on `(0,1,0)` alone and ~216 on all six.
- **Best known parameter set: `J1 = 0.15`, `sigma_J = 0.50`, `gzz = 3.50`,
  `sigma_gzz = 0.80`.** `chi2_red = 24.1` across all six cuts against `110.4` for the by-eye
  set, per-cut 10.6 to 50.5. **Both disorder widths land exactly on the by-eye values** while
  `J1` falls 40% and `gzz` 8%. So the by-eye disorder was right, and what needed correcting
  were the two parameters inherited from the analytical model -- which cannot represent mode
  mixing between exchange environments, whereas Sunny/KPM can. Reached by a coordinate
  descent that was still improving when stopped, so it is the best KNOWN point, not a
  converged optimum. One full six-cut evaluation at 81 q with 4 realizations is 456 s.
- **Away from the zone centre the dispersion runs DOWNWARD from the Zeeman energy**, so
  lowering `J1` *raises* the K and M peaks toward `gzz*mu_B*B`, while lowering `gzz` lowers
  everything. At K and 9 T: by-eye 1.98 - 1.13 = 0.86 meV, fitted 1.82 - 0.68 = 1.15 meV,
  both matching the observed peak positions. That is why the by-eye set could not fit both
  ends at once -- both parameters were too large, pushing Gamma too HIGH and K/M too LOW
  simultaneously. An earlier note here read that pattern as "the bandwidth is too small",
  which was exactly backwards.
- **RETRACTED: "the exchange sector fails at K and M".** An earlier version of this file
  recorded `chi2` of 164-285 at K/M and concluded that isotropic Heisenberg cannot describe
  the dispersive cuts, with XXZ the likely remedy. That was an artefact of a scan grid which
  fixed `gzz = 3.30` and floored `J1` at 0.20, while the good region is `gzz = 3.50` with
  `J1 = 0.15` -- outside the grid on both axes. At the right parameters K/M sit at 12-50 and
  the Gamma-to-K/M gap is 2-3x, not 25x. **Do not revive the XXZ argument from a chi2
  surface without first checking that the grid contains the optimum.**
- **The M(H) linear term is not Van Vleck, and this is now settled three ways.** Our crystal
  field gives `chi_VV^zz = 0.0171 +- 0.0007 uB/T`. The **Bag et al. Supplement independently
  fits 1.39(2)e-2** for H ∥ c at 2.5 K — agreeing with our crystal field to ~25% by a wholly
  different route. Both are **~13x below the 0.19 uB/T** our M(H) fit wanted, and forcing that
  term to zero degraded rms from 0.0107 to 0.0652, so the fit leaned on it hard. The excess is
  therefore neither Van Vleck nor a bad Van Vleck estimate: it is **our model's approach to
  saturation being the wrong shape**, with the linear term patching it. Consistent with the
  neutron optimum overshooting `B_sat` to ~3.3 T where M(H) wants ~4.0 T.
- **We AGREE with the published work on the g factor, so the argument is not about g.** The Bag
  et al. Supplement fits **g_par = 3.436(4)** at 2.5 K, which sits inside our clean neutron
  determination (`gzz` ~ 3.40, Gamma surface 3.30-3.45). The disagreement is about the
  **exchange** and about **disorder**, not about g. State it that way — it narrows the argument
  to ground this work is strong on.
- **Our sample reproduces theirs to 4.2%, identically in both orientations.** At 14 T and 2.5 K
  our DynaCool VSM gives 1.84 (∥) and 1.59 (⟂) uB/Yb against the Supplement's ~1.92 and ~1.66 —
  ratio 0.958 in *both*. A uniform offset across orientations is a normalization systematic
  (mass, holder, demagnetization), not a sample difference. So the sample behind the erroneous
  claim and ours show the same observables.
- **The g tensor is nearly ISOTROPIC, which argues against putting all anisotropy in g.**
  The Supplement fits **g_perp = 3.037(5)**, so `g_perp/g_par = 0.884`. Our own 2.5 K data then
  point two ways: at 14 T `M_par/M_perp = 1.157` against `g_par/g_perp = 1.131` (consistent, so
  high-field anisotropy *is* g anisotropy), but at 1 T `M_par/M_perp = 1.50` against the ~1.28
  that `g^2` predicts (**not** consistent). Extra anisotropy appearing where correlations matter
  and vanishing where the polarized state dominates is the signature of **exchange** anisotropy.
  So some XXZ may be genuinely needed — on top of the disorder model, not instead of it.
  Caveat: at 1 T and 2.5 K this is not linear response (`g mu_B B` ~ 0.20 meV vs
  `kT` = 0.215 meV), so it is suggestive rather than settled. **This is a finding against the
  modelling choice in this repo and is recorded deliberately.**
- **Where the neutron fit's chi2 comes from, per energy.** The 1.8-2.4 meV band is 24% of the
  fit window's width but carries 45% of `chi2` on average. It must be read against where each
  mode sits: `(0,1,0)` 9 T is 72% but its mode is *in* the band and that cut carries no 2.08 meV
  background, so that share is the magnon and is desirable. The red flags are the 9 T dispersive
  cuts at **35-41% with their modes at ~1.05-1.10 meV, outside the band** — that share is
  instrumental artefact, and it is the measured mechanism behind the gzz pull to 3.70.
  `(0,1,0)` 14 T is only 16% because its mode has essentially left the window, which is the
  window-exit artefact seen from another angle.

## Protocol for M(H) work

T = 0 `minimize_energy!` from a field-polarized start with adiabatic field
continuation, 12x12x1, 8-16 **fixed** disorder realizations (common random numbers)
so the objective is deterministic in the parameters. Validate any optimum at
36x36x1 with different seeds. Threading here is over REALIZATIONS rather than over q, so
the KPM thread-count warning below does not apply directly. But `-t auto` is still the wrong
choice on a many-core box, because threads beyond the realization count do nothing at all:
use `julia -t N` with N equal to the number of realizations (8-16). The realization path's
own scaling knee has not been measured.

## Gotchas

- **`@printf`/`@sprintf` need a LITERAL format string.** `@printf("a" * "b", x)` fails with
  "First argument ... must be a format string", and it fails at MACRO EXPANSION during load,
  so it survives `Meta.parseall` and dies only when the script runs -- after the compute. This
  has cost three separate cycles here, one of them 950 s of KPM. Wrap long formats by putting
  the ARGUMENTS on continuation lines, never by concatenating the format itself.

- **Do NOT run KPM under `julia -t auto` on a many-core box.** Intra-process q-threading
  loses efficiency fast, and past a knee it goes *negative*. Measured intra-process
  efficiency (one 81-q spectrum, 36x36x1):

  | chunks/threads | 1 | 2 | 4 | 8 | 16 | 32 | 128 |
  |---|---|---|---|---|---|---|---|
  | Windows, 32 cores | 100% | 92% | 73% | 36% | 22% | 10% | - |
  | DGX, 2x EPYC 7742 | - | - | ~100% | 90% | 73% | - | **5%** |

  On the DGX, **one process with 128 threads is the worst configuration measured** -- 21.96 s
  per spectrum against 11.25 s for the *same code* on 16 threads. The fix is to fan out over
  **processes**, threading only as far as it stays efficient: 32 processes x 4 threads gives
  0.859 spec/s against 0.046 for `-t auto`, an **18.7x** difference. On the Windows box the
  equivalent is ~8 processes x 4 threads, worth ~7x over the 3.4x that 32 threads saturates at.
  The ceiling is process-internal (synchronisation or per-chunk overhead in the q-loop), NOT
  memory bandwidth: per-process throughput stays flat as processes are added, 8x16 and 16x8
  differ by 17% at identical total threads, and NUMA pinning made it *worse*. An earlier claim
  in this file that KPM is "memory-bandwidth bound" was wrong.
  Corollaries: keep total threads at or below the PHYSICAL core count (using SMT siblings
  costs throughput); do not `taskset`-pin (it created stragglers on the nodes holding the
  shared Julia sysimage page cache); use dynamic work assignment, since at 8x16 the batch span
  was 105.6 s against a 70 s mean and one straggler paces everything.
- **`relax_attempts` must stay at 1.** The escalation loop in `sv_kpm_context` breaks on
  `minimize_energy!`'s convergence *flag*, not on a KPM rejection, and at 9 T the flag never
  trips -- so anything above 1 escalates maxiters 1000 -> 4000 -> 16000 unconditionally,
  chasing exactly what this file says to ignore. Measured: **21.8x the relaxation cost for
  5.8e-9 in E/site**, and the spectrum changes by at most 0.12%, ~100x below the realization
  floor. `sv_neutron_curves` defaulted to 3 until 519cf2f and silently cost ~1.55x on every
  neutron evaluation. Raise it only if KPM actually rejects states -- that is a
  *regularization* problem, not a maxiters problem.
- `sv_interp1` returns NaN outside the measured range. The M(H) data span
  0.0225-6.975 T; the objective configs use 0.2-6.8 T for margin, which discards
  usable low-field data (see open threads).
- The KPM kernel FWHM is bounded above by the instrument: the CNCS table *falls*
  with energy transfer (0.155 meV at 0.5 meV to 0.055 at 4 meV). A kernel wider
  than the target cannot be corrected by deposition. Use
  `sv_resolution_kernel_headroom` before a run.
- `[kpm].tol = 0.05` carries ~10% rms error against tol = 0.005. Acceptable for
  looking at a lineshape, not for fitting; 0.01 costs only 1.7x.
- The MC Q-sampling mode at `n_events = 5000` is ~36 min per 6-cut comparison on
  Windows CPU but only ~4 min per cut on one A100, so it is no longer cost-blocked.
  Whether the extra events buy anything over the 81-point grid is still untested.
- `results/` is gitignored, so figures do not survive a clone. Regenerate them.

## Open threads

1. **SETTLED - q sampling: use 81 q (3x3 measured x 3x3 Gauss-Hermite resolution).** The
   deterministic grid, not MC. This REVERSES an earlier recommendation of 225 q here, which was
   convergence for its own sake: 225 q is better converged (0.058% against 0.64% versus a 625-q
   reference) but the two differ by only 0.286 in `chi2_red` at n = 8, well below the ~1.0
   realization scatter and far below the background systematic, while costing 1.76x more. Extra
   *resolution* nodes buy nothing at all -- three Gauss-Hermite nodes are already exact -- so the
   *measured* axis is the only one that matters. MC is not preferred: at 625 events it still
   carries ~0.5% sampling error, and while its samples are frozen across evaluations (hence
   deterministic) they are irregular in parameter space. All of this sits far below the systematic
   floor described above, which is the general lesson: check where the floor is before buying
   numerical precision.
2. **MOSTLY DONE - the neutron objective.** `sv_neutron_objective` exists, with realization
   averaging under common random numbers, q-threading, a profiled-out intensity scale,
   failure sentinels and per-cut reporting. What is still missing is an OPTIMIZER: the scans
   use grid search plus coordinate descent driven from a script. Also outstanding is the
   realization scatter expressed in `chi2_red` units rather than spectrum units, which is the
   number a fit actually consumes.
3. **Co-optimization** of the spectra and M(H), which is the point of all of the
   above: M(H) pins B_sat and sigma_gzz, the spectra must pin sigma_J.
4. Does the neutron model still need the phenomenological flat component (r2)
   at high disorder? M(H) no longer does.
5. Where does the excess linear M(H) term come from, if not Van Vleck? Suspect the
   normalization of the digitized data, an impurity, or model error.
6. **Co-fit the background with the model refinement.** Noted deliberately, NOT pursued yet.
   The background is currently constructed first and frozen, then the model is fitted to the
   corrected data -- so background error propagates into parameters with no route back. Fitting
   a parameterised background jointly with the exchange and g parameters would let the data
   constrain both, and would produce a proper covariance between them rather than a one-way
   systematic. The obvious hazard is that a flexible background can absorb real signal, so it
   would need tight physical priors (the A1-A4 assumptions in
   `scripts/background_stage1_ei332.jl`) and a demonstration that it does not simply eat the
   magnon. Worth doing after the per-point variance route is exhausted, not before.
7. **Non-Gaussian disorder distributions.** If quantitative agreement keeps failing, the
   next thing to try is a non-Gaussian distribution of `gzz` and/or `J`, NOT more XXZ.
   There is a physical reason to expect this: Zn/Ga site mixing gives each magnetic site a
   *discrete* count of Ga vs Zn neighbours, so the natural distribution is multinomial over
   local environments -- multi-modal, not a smooth Gaussian. A Gaussian is the convenient
   parameterisation, not the physically motivated one. Deliberately deferred, not dismissed.
8. Zero field. The published claim is a zero-field statement and everything here is
   field-polarized. Hard: no polarized reference state, LSWT invalid, and classical
   statistics is a poor approximation at the 0.07 K neutron temperature.

## Two-machine workflow — read this before touching git

This repo is worked on from **two machines that share state only through GitHub**:

- a **Windows desktop** (paths under `C:\Users\vdp\ORNL Dropbox\...`), and
- **`neutrons-dgx01.ornl.gov`**, a *shared* Linux box with A100s, at
  `~/repos/YbZn2GaO5`, driven over VS Code Remote-SSH.

Remote is `github.com/dpaj/YbZn2GaO5`; work goes straight to `main`. There is no
direct link between the two sessions, and **two separate Claude Code sessions run
on the two machines and cannot see each other.** The human is the only channel
between them, and this file is the only shared memory — per-machine Claude memory
directories do not travel.

Rules that follow from that:

- **`git fetch` before reporting or reasoning about repo state.** `git status`
  saying "up to date with 'origin/main'" only compares HEAD to the *locally cached*
  remote ref. A box that has not fetched will report "up to date" while being many
  commits behind. This has already caused confusion once.
- **Commit and push before switching machines or ending a session.** Do not leave
  work stranded on one side.
- **Never assume the local working copy is current** if the other machine may have
  been touched. Ask.
- **Prefer `git pull --ff-only`.** It refuses rather than silently creating a merge
  commit when the two sides have diverged, which is what you want when you cannot
  see the other machine.
- When both sides have work: **pull first, then commit.** Committing first makes the
  box simultaneously ahead and behind for no reason.

### Files that need care across machines

- **Set your git identity on every machine**, or git invents one from
  username@hostname (which is how `vdp@neutrons-dgx01.ornl.gov` reached history):

      git config --global user.name  "Daniel Pajerowski"
      git config --global user.email "daniel@pajerowski.com"

  That address is what most existing commits use and what links to the `dpaj`
  GitHub account. `.mailmap` canonicalizes the historical variants for git tooling,
  but GitHub does not honour `.mailmap` for account linking, so setting this up
  front is the only real fix.
- `CLAUDE.md` and `.claude/settings.json` are **shared project config** — push and
  pull them like code. `.claude/settings.local.json` is machine-specific and is
  ignored by this repo's `.gitignore` (deliberately not by a machine-local global
  ignore, which would not travel).
- **`Manifest.toml`, `envs/sunny-main/Manifest.toml`, `envs/sunny-kpm-gpu/Manifest.toml`
  are tracked and will churn across OSes.** They carry 167-185 `_jll` deps, and the
  CUDA artifacts in the GPU env resolve differently on Linux. `julia_version` is
  pinned at 1.12.3, so a different Julia guarantees a diff. After running
  `Pkg.instantiate()` on a second machine, **do not commit the churn** — use
  `git checkout -- Manifest.toml` — or it ping-pongs forever.
- There is **no `.gitattributes`**, and 111 tracked files carry CRLF in the index.
  Renormalizing is a one-time noisy commit touching all of them, so it must be done
  while *both* machines are clean and pulled immediately on the other side.
- Three `scripts/legacy/*.jl` files have hardcoded `C:\Users\vdp\...` paths and will
  fail on Linux. They are still in use, so the fix is to read an env var with the
  Windows path as fallback, following the pattern already at
  `plot_yzgo_2d_data_vs_model_legacy.jl:5910`
  (`get(ENV, "YZGO_DATA_DIR", raw"C:\...")`). Note the data they want lives in
  `YZGO/CNCS_data/`, **outside this repo**, so it must be copied to the DGX
  separately.
- `../references/` (13 MB of published PDFs) is outside the repo and must stay
  untracked — it is a shared box, and they are copyrighted.

### Which machine for what

Windows box: 32 cores, RTX A2000 (6 GB, FP64 at 1/32 rate). CPU q-threading
saturates near 3x here. DGX: 8x A100-SXM4-40GB, of which **GPUs 4-7 are typically
held by vLLM** (~0.7 GiB free) and 0-3 are usable; CPU q-threading reaches 16.6x.

Measured on the DGX at 36x36x1, 81 q, tol 0.05, kernel 0.05 meV:

**Caveat on every figure in this table and the 3x/16.6x threading numbers above:** they were
measured with `regularization` at Sunny's 1e-8 default, not the production 1e-5, because the
benchmark scripts omitted the keyword (`benchmark_kpm_1d_scaling.jl`,
`benchmark_sunny_kpm_core.jl`, `benchmark_sunny_kpm_cpu_gpu_compare.jl`). Regularization changes
the spectral range and hence the Chebyshev moment count, so these are correct measurements of a
*non-production* configuration. Treat them as relative comparisons between paths, not as absolute
production costs. The production per-evaluation cost, measured properly, is **466 s for a six-cut
objective at 81 q with 8 realizations on 32 DGX threads**.

| path | s per q |
|---|---|
| CPU serial | 0.774 |
| CPU, 81 chunks | 0.0495 |
| GPU Float64 | 0.0653 |
| GPU Float64, 441-q batch | 0.0488 |
| GPU Float32 | 0.0532 |
| 4 GPUs concurrent | 3.91x aggregate |

So one A100 in Float64 roughly ties well-threaded DGX CPU; **the win is running
several GPUs at once**. Device memory is ~0.85 MiB per q, so the 5000-event MC mode
needs ~4.8 GiB -- impossible on the 6 GB A2000, trivial on a 40 GB A100 (~45,000 q
per card). Heavy KPM and the MC Q-sampling mode belong on the DGX; M(H) work is
cheap enough anywhere. Note the DGX CPU serial figure is not comparable to the
Windows box 0.216 s/q: it was measured with the kpm-gpu fork of Sunny 0.9.0 and on
slower per-core hardware.

## Housekeeping

Run `julia --project=. test/runtests.jl` before committing; it catches convention
regressions. On a machine that has just pulled, also check `julia --version`
against the manifests' 1.12.3 and run `Pkg.instantiate()` first.

**Do not quote a fixed assertion count.** The total is not a constant: the config
testset asserts once per file in `configs/`, and the KPM-threading testset
self-skips without threads (contributing two assertions with `-t auto`, zero
without). So the number moves whenever a config is added or an assertion is written,
and it differs between machines and between threaded and unthreaded runs — observed
values include 53, 54, 55, 56 and 57, all correct. **Judge it by "all green, 0
failures".** Runtime is machine-dependent too (about 20 s on the Windows box, 34 s
on the DGX), so it is not a regression signal either.
