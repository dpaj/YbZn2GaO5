# Large-cell classical M(H,T): minimal single-disordered-phase model

Workflow added by:

```text
scripts/sunny_largecell_mvh_classical.jl
configs/sunny_largecell_mvh_classical_controls.toml
```

Run with:

```powershell
julia --project=. scripts/sunny_largecell_mvh_classical.jl
```

A cheaper or alternative control file can be substituted without editing the
production config, either as a positional argument or via `SUNNY_MVH_CONTROLS`:

```powershell
julia --project=. scripts/sunny_largecell_mvh_classical.jl configs/my_variant.toml
```

## What this calculation is

Magnetization versus field at finite temperature for a **single disordered
phase**, evaluated on a tunable Sunny supercell in the same effective-Hamiltonian
scheme as the KPM LSWT neutron calculation:

- the same effective P1 one-site triangular net (`sv_effective_triangle_crystal`),
- the same explicit J1 and J2 bond shells,
- the same fractional per-bond exchange disorder `J*(1 + sigma_J*randn)`,
- the same per-site `gzz` disorder,
- the same field direction (H parallel to c) and the same base seed,
- `dims = [3,3,1]` seed cell enlarged by `repeat_periodically`, so
  `cell_sizes = [[36,36,1]]` reproduces the KPM neutron supercell exactly.

Two things are deliberately different from the analytical co-fit.

**One phase, not two.** The separate non-dispersive / "flat" component
(`gzz2`, `sigma_gzz2`, weight `r2`) is not included. The premise under test is
that thermal population of the low-energy modes that disorder generates inside a
single *coupled* Hamiltonian does the work that the analytical model does with a
second, independent, non-dispersive phase. This is the M(H) analogue of the
`README.md` distinction:

> Analytical sigma_J: distribution of clean dispersions from different effective regions
>
> Sunny bond sigma_J: one disordered Hamiltonian with many exchange values coupled together

The reported model is therefore

```text
M_total(B) = A_M * [ M_disp(B,T) + chi_vv * B ]
```

with no `r2` mixing and no `1/(1+r2)` normalization. Note this differs by a
factor `(1+r2) = 1.158` from the normalization used by the other Sunny
magnetization paths, so `A_M` is not directly comparable between them. The free
least-squares `A_M` is reported alongside the fixed one for this reason.

Van Vleck is a single-ion term on the same Yb rather than a second phase, so it is
kept in the minimal model. Note `A_M` multiplies it too, so the contribution to
the observable is `A_M * chi_vv * B`, which at `A_M ~ 0.52` is 0.29 uB at 7 T
against a measured 1.12 uB — about a quarter of the signal, not the 0.55 uB the
unscaled slope alone suggests.

**But `chi_vv_muB_per_T` was fitted inside the analytical model and carries that
model's assumptions.** Holding it fixed here can manufacture a residual that looks
like missing physics, so it is fitted jointly with `A_M` by nonnegative least
squares and reported as `chi_vv_joint_fit_muB_per_T`. In the second production run
below it fits to **exactly zero** in every case, which halves the residual. Do not
carry the analytical value into this context without checking it.
`include_chi_vv = false` drops the fixed-chi_vv convention entirely.

**Finite temperature by classical sampling, not T = 0 minimization.** The
existing `sv_sweep_largecell_component` path is `minimize_energy!` only.

## Samplers

Sunny 0.9.1 has exactly two integrators and only one of them is a thermostat:

| Integrator | Role here |
|---|---|
| `Langevin(dt; damping, kT)` | stochastic Heun; samples the classical Boltzmann distribution. The only thermal sampler. |
| `ImplicitMidpoint(dt; tol)` | symplectic and energy conserving. **Errors if `kT != 0`** (`Integrators.jl:139`), so it cannot thermalize. |

The `samplers` list therefore selects among:

- `minimize_energy` — T = 0 reference, i.e. the existing large-cell scheme.
- `langevin` — relax at each field, then heat to `kT` and time-average.
- `langevin_then_midpoint` — relax, thermalize with Langevin, then time-average
  under conservative dynamics. This removes the thermostat's discretization from
  the sampled observable and is an independent cross-check rather than a
  different physical ensemble.

`relax_before_thermalize = true` is not cosmetic. Thermalizing directly from the
previous field's thermal configuration does not converge in any affordable
number of steps: a test run without it returned `M(7 T) ~ 0.02 uB` instead of
`~1.73 uB`, and the microcanonical variant went negative. Each field is relaxed
with `minimize_energy!` first, then heated.

## Known systematic: classical over-counting of high-energy modes

Classical statistics assigns every magnon mode an occupancy `kT/eps`, where
quantum statistics gives `1/(exp(eps/kT) - 1)`. These agree when `eps << kT`,
which is exactly the low-energy disorder-generated regime of interest. They
disagree badly when `eps >> kT`, where the classical result over-counts.

The validation panel measures this directly. Above saturation the field gaps the
spectrum, so the quantum answer is the fully saturated moment to within about
`1e-8`, and any Langevin deficit there is pure classical error. Measured on this
Hamiltonian at 0.42 K, 36x36x1, 9 T:

| Quantity | Value |
|---|---|
| minimum magnon energy | 0.421 meV |
| `kT` at 0.42 K | 0.0362 meV |
| `sum n_B` over 1296 modes | 1.4e-5 |
| LSWT + Bose moment reduction | 4e-8 uB/site |
| classical Langevin deficit | ~0.11 uB/site (~6%) |

The deficit is field-dependent, so it distorts the **shape** of M(H) and not
merely its scale. `validation_fields_T` and `validation_lswt` control this panel;
it costs one dense `2N x 2N` diagonalization per field (about 10 s at N = 1296).

Sunny's `set_spin_rescaling_for_static_sum_rule!` (`|S|^2 = s(s+1)`) is exposed
as `spin_rescaling_for_static_sum_rule` but **defaults to false**, because
Sunny's own documentation states that this rescaling is a high-temperature fix
and is explicitly not appropriate at low temperature.

## Why LSWT + Bose is not used for the 0-7 T curve

Finite-temperature LSWT would be the exact finite-T version of the KPM
approximation scheme, and `SpinWaveTheory` does accept the inhomogeneous 1296-site
system (1296 bands, about 10 s per diagonalization). But the model is strictly
two-dimensional Heisenberg with no single-ion anisotropy, so below saturation the
spectrum is gapless and the magnon occupation sum diverges. Measured minimum
magnon energies on the 36x36x1 cell: 4.8e-5 meV at 4.5 T and 8.5e-5 meV at 1 T,
giving `sum n_B` of 839 and 465 respectively, i.e. an unphysical moment reduction
larger than the total moment. The measured M(H) window is entirely at or below
saturation, so the LSWT route is retained only for the above-saturation
validation panel.

## Cell size is a measurement, not a knob

Because the model is strictly 2D with no interlayer exchange, below saturation
there is a U(1) Goldstone mode, and Mermin-Wagner implies the classical thermal
moment reduction should grow like `log(L)` rather than converge. Above saturation
the field gaps the spectrum and convergence is expected.

Changing the cell size also changes which random numbers land on which bond, so a
cell-size scan is only interpretable when averaged over disorder realizations.
`sv_apply_disorder!` gained an optional `realization` keyword for this;
`realization = 0` reproduces the original single-realization seed bit-for-bit, so
existing neutron and KPM results are unchanged.

Report `M` versus `L` rather than selecting one size. If the low-field curve
drifts with `L`, the remedy is a small interlayer exchange, not a larger box.

## Convergence checks

Equilibrium averages must be independent of `damping`. If they are not, `dt` is
too large. `damping_scan` and `dt_scan` run this at `convergence_field_T` and
write `..._convergence.csv`. On this Hamiltonian `suggest_timestep` wants
`dt ~ 0.078` at `tol = 1e-2`, essentially independent of damping; the config
default `dt = 0.02` is deliberately conservative. `dt_mode = "suggest"` calls
`Sunny.suggest_timestep_aux` and applies `dt_safety`.

## Outputs

```text
results/feature_tables/sunny_validation/largecell_mvh_classical/
    sunny_largecell_mvh_classical.csv             per-point, per (cell, realization, sampler)
    sunny_largecell_mvh_classical_manifest.csv    one row per sweep, all knobs
    sunny_largecell_mvh_classical_profile.csv     per-stage timings
    sunny_largecell_mvh_classical_validation.csv  above-saturation calibration
    sunny_largecell_mvh_classical_convergence.csv dt / damping scan
results/figures/sunny_validation/largecell_mvh_classical/
    sunny_largecell_mvh_classical.png
```

Every row records the Sunny version, cell size, seed, realization index,
sampler, `dt`, `damping`, step counts, temperature, and both the fixed and free
`A_M`, per the provenance requirement in
`docs/companion/05_sunny_validation.md`.

`M_sat_realization_uB` is `mean(g_i) * S` for that specific realization, not the
nominal `gzz * S`. Per-site `gzz` is disordered, so a small cell can legitimately
exceed the nominal saturation value by `O(sigma_gzz/sqrt(N))`; deficits are
measured against the realization value.

## Absolute scale

The primary curve uses `magnetization_global_scale` from
`configs/best_fit_parameters.toml`. The least-squares `A_M` is reported next to
it rather than replacing it, so that the absolute-moment discrepancy stays
visible instead of being absorbed into a scale factor. This matters: the model
saturates near 5 T at `gzz*S = 1.892 uB/site`, whereas the measurement is at
`1.117 uB/Yb` at 6.97 T and still rising. `A_M = 0.419` is absorbing a factor of
roughly 2.4 in absolute moment, which is plausibly the same issue as the
"weirdness for the phase fractions" note in `README.md`.

## First production run

Config as committed: 0-7 T in 36 steps, 0.42 K, cell sizes 12/24/36, three
disorder realizations, all three samplers. 2508 s of profiled wall clock.

Cost by sweep set (3 realizations x 36 fields):

| Cell | `minimize_energy` | `langevin` | `langevin_then_midpoint` |
|---|---|---|---|
| 12x12x1 | 11 s | 41 s | 117 s |
| 24x24x1 | 39 s | 211 s | 370 s |
| 36x36x1 | 178 s | 541 s | 979 s |

**Sampler cross-check passes.** Equilibrium `M` is independent of `damping` to
0.0005 uB at `dt = 0.02` (0.05/0.1/0.2 give 1.78164/1.78118/1.78143), well inside
the 0.0025 uB sampling error. `langevin` and `langevin_then_midpoint` agree to
0.002 uB, so the thermostat is not biasing the observable. `dt = 0.01` with
`damping = 0.05` is the one mild outlier, consistent with slower decorrelation at
fixed sample count rather than a discretization problem.

**The large box is not needed for M(H).** The free `A_M` is 0.5225, 0.5248, 0.5225
for 12/24/36 with `langevin`, and 0.4727, 0.4741, 0.4720 with `minimize_energy` —
flat to under 0.5%. The predicted `log(L)` Goldstone drift did **not** appear.
This is expected in hindsight: `M_z` along the field is not the symmetry-broken
order parameter, so transverse Goldstone fluctuations do not degrade it. So
`12x12x1` is sufficient for magnetization, and a 12x12 Langevin sweep at 36 fields
costs about 14 s per realization, which makes fitting against M(H) affordable.

**Classical over-counting, measured.** 0.111 uB (5.88%) at 9 T and 0.058 uB
(3.06%) at 14 T, against an LSWT + Bose quantum answer of 4e-8 and 4.5e-18 uB.
Field-dependent and shrinking as the gap opens, as expected.

**Finite temperature barely matters at 0.42 K, and does not supply the missing
low-field curvature.** Maximum residual against experiment at free `A_M` is
0.062 uB for T = 0 versus 0.069 uB for classical 0.42 K — i.e. thermal sampling
changes M(H) by less than the model-versus-experiment systematic, and slightly
*worsens* it. Both leave the same S-shaped residual on the 36x36x1 cell:

| B (T) | experiment | model (free A_M) | residual |
|---|---|---|---|
| 0.6 | 0.2137 | 0.1452 | -0.069 |
| 1.0 | 0.3139 | 0.2487 | -0.065 |
| 2.0 | 0.5266 | 0.4731 | -0.054 |
| 3.0 | 0.7249 | 0.6883 | -0.037 |
| 4.6 | 0.9734 | 0.9723 | -0.001 |
| 5.8 | 1.0685 | 1.1019 | +0.033 |
| 6.8 | 1.1111 | 1.1756 | +0.065 |

The model is too straight: it undershoots below the 4.7 T crossover and
overshoots above it. Experiment has both more initial susceptibility and a
stronger rollover. That is the signature of a missing fast-saturating component,
which is what the analytical model's `tanh` non-dispersive phase supplies.

**Why the minimal model cannot produce it as configured.** The bond disorder is
multiplicative, `J*(1 + sigma_J*randn)`, so at `sigma_J = 0.24` it rescales bonds
but never decouples a site. Genuinely free or near-free moments — the states whose
thermal population would give low-field curvature — require either a much larger
`sigma_J` (cf. the 2-3x multipliers already being explored on the neutron side in
`configs/sunny_kpm_1d_disp_grid_2sigmaJ_controls.toml`) or a different disorder
channel such as bond removal / site dilution. Note also that `sigma_gzz` disorder
is clipped at `max(0, ...)`, so low-`g` sites become *invisible* rather than free.

**Absolute moment.** Free `A_M` of 0.47-0.52 means the model moment is close to
twice the measurement. Combined with the model saturating near 5 T while the data
are still rising at 6.97 T, this points at `gzz` being too large for the
magnetization, which is plausible given that the co-fit that produced
`gzz = 3.785` has neutron reduced chi-squared of 16.6.

## Second production run: by-eye KPM parameters, sigma_J = 0.5

`[run.param_overrides]` set to the by-eye Sunny-KPM neutron parameters
(`J1 = 0.25`, `J2 = 0.01`, `sigma_J = 0.5`, `gzz = 3.8`, `sigma_gzz = 0.8`), with
a `clean` variant (both disorder widths zeroed) overplotted and `chi_vv` fitted
jointly with `A_M` instead of held at the analytical value. 1041 s.

**Disorder produces both of the effects it was expected to, and they are T = 0
effects.** Raw moment per site, 36x36x1, `minimize_energy`:

| B (T) | disordered | clean |
|---|---|---|
| 0.2 | 0.166 | 0.074 |
| 1.0 | 0.507 | 0.372 |
| 5.0 | 1.617 | 1.857 |
| 5.4 | 1.689 | **1.9000** |
| 7.0 | 1.856 | **1.9000** |

The clean system pins at exactly `gzz*S = 1.9000` from 5.4 T upward — the hard
kink. The disordered system rounds it off and never saturates within the window.
And the initial slope more than doubles at 0.2 T. Both are the behaviours the
old implementation showed.

**Fit quality, 36x36x1** (RMS and max |residual| against experiment, uB/Yb):

| variant | sampler | fixed analytical chi_vv | chi_vv fitted jointly |
|---|---|---|---|
| disordered | minimize_energy | 0.0354 / 0.0779 | **0.0182 / 0.0352** |
| disordered | langevin | 0.0377 / 0.0822 | **0.0176 / 0.0365** |
| clean | minimize_energy | 0.0717 / 0.1012 | 0.0665 / 0.0970 |
| clean | langevin | 0.0608 / 0.0915 | 0.0496 / 0.0821 |

Two independent factor-of-two improvements: turning the disorder up, and letting
`chi_vv` float. RMS 0.018 uB is about 1.6% of the 1.12 uB signal, against
max 0.062-0.069 uB for the canonical `sigma_J = 0.24` set. Disorder is essential
— it improves RMS by a factor of about 3.7 over clean.

**chi_vv fits to exactly zero in all eight (cell, sampler, variant)
combinations.** The nonnegative joint fit could have selected the analytical
`chi_vv = 0.0783` — it is inside the feasible set — and instead drives it to the
boundary while raising `A_M` from 0.540 to 0.685. So the analytical Van Vleck
slope was absorbing high-field curvature that the coupled disordered Hamiltonian
now supplies directly. It should not be carried over into the Sunny M(H) context
unexamined.

**Temperature still does not matter; disorder does.** T = 0 and classical 0.42 K
give RMS 0.0182 versus 0.0176, i.e. indistinguishable. At low field the T = 0
disordered curve actually has *more* initial slope than the classical finite-T one
(0.166 versus 0.107 uB at 0.2 T), so classical thermal fluctuations slightly
*reduce* the low-field moment. The extra initial susceptibility is a T = 0
disorder effect — weakly coupled spins that turn over in a small field — not
thermal population of soft modes.

**Cell size remains irrelevant:** `A_M_joint` is 0.6847 at 12x12x1 versus 0.6845
at 36x36x1.

**Absolute moment still off by about 1.5x.** `A_M = 0.685`, improved from about
0.42-0.52 but still far from 1. That corresponds to an effective
`gzz ~ 3.8 * 0.685 = 2.6` for magnetization rather than 3.8.

**Residual shape that survives:** mildly negative (about -0.03 uB) below 2.5 T
and positive (+0.037 uB) above 5.5 T. The S-shape is much reduced but not gone.

**Validation note.** With `sigma_gzz = 0.8` the disordered system is *not*
saturated at 9-14 T, because low-`g` sites need many times more field, so the
LSWT reference is invalid there — during testing it reported a nonsense 27.9 uB
depletion. The calibration therefore runs on the `clean` system
(`validation_variant`), where 9 T gives a gap of 23.6 kT and LSWT + Bose confirms
a quantum depletion of 1.8e-11 uB. Classical over-counting there is 0.116 uB
(6.11%) at 9 T and 0.059 uB (3.12%) at 14 T. The script now flags any validation
point whose gap is not comfortably above kT rather than printing a bad number.

## Open items

- The measurement temperature is unreconciled: 0.42 K in the control TOMLs,
  0.4 K as a legacy co-fit keyword default, and no recorded provenance for
  either. `docs/companion/04_magnetization_model.md` still lists temperature as
  an unfilled item.
- The M(H) window (0-7 T) does not overlap the neutron fields (9 T, 14 T).
- Whether a second phase is needed at all is the question this workflow exists
  to answer; if the minimal model fails, the flat component can be reintroduced
  through `sv_build_magnetization_comparison`.
- Bose factors are still absent from the KPM LSWT neutron calculation
  (`intensities` is called without `kT`). The effect is expected to be small at
  0.07 K but should eventually be included.
