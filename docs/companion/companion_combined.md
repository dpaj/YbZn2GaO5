# YbZn2GaO5 analysis companion document


This combined Markdown file concatenates the modular companion docs. For editing, prefer the individual files in `docs/companion/`.




## Companion documentation overview

This directory is a starter companion-document layer for the `dpaj/YbZn2GaO5` GitHub repository. It is written as modular Markdown so that each section can mature independently while still forming a coherent reproducibility guide.

## Recommended structure

- `01_experiment_and_run_inventory.md` — sample photographs, neutron experiment inventory, run-condition summary, and data provenance.
- `02_neutron_reduction_backgrounds.md` — inelastic neutron cuts, histogramming, background subtraction, and 3.32/4.65 meV comparison.
- `03_analytical_field_polarized_neutron_model.md` — analytical field-polarized neutron model, disorder averaging, and experimental observable generation.
- `04_magnetization_model.md` — mean-field/disorder magnetization model and co-fit implementation notes.
- `05_sunny_validation.md` — Sunny.jl validation layer, current approximations, and comparison to the analytical model.
- `06_reproducibility_and_provenance.md` — script/config/figure mapping, run instructions, and release checklist.

## Current maturity map

| Topic | Status | Recommended next action |
|---|---|---|
| Neutron data collection and run histogram inventory | mature enough to document | Add final histogram figures and link raw/reduced data products. |
| Background subtraction and 3.32/4.65 meV comparison | mature enough to document | Add final comparison figure filenames and state which background mode is canonical. |
| Analytical field-polarized neutron model at 9 T and 14 T | mature enough to document | Add derivation equations checked against implementation. |
| Analytical model to neutron observables | mature enough to document | Document bin-center versus full histogram treatment explicitly. |
| Magnetization model and implementation | mature enough to document | Add units-conversion details and final M(H) component plots. |
| Sunny validation | active / partially mature | Keep separate from the primary analysis until approximation choices stabilize. |

## Tables generated from the attached run CSV

- `tables/neutron_run_summary_core_analysis.csv`
- `tables/neutron_run_summary_by_condition.csv`
- `tables/neutron_run_inventory_parsed.csv`

The parsed run CSV covers 4482 runs for sample `YbZn2GaO5`, run numbers 663980–668461, with total listed proton charge 401.46 C, from 2025-10-17 18:11:27 to 2025-10-27 08:09:13. The sample environment listed in the run table is `14T Vertical Field Magnet;Oxford Dil Fridge insert`.




## 1. Experiment and run inventory

## Purpose

This section documents the experimental data products that underlie the YbZn2GaO5 neutron and magnetization analysis. The goal is to make it clear which measurements were performed, which subset is used in the current modeling workflow, and how the run metadata maps onto reduced data files and analysis scripts.

## Sample and crystal photographs

Add crystal/sample photographs here.

Suggested files:

| Placeholder | Suggested path | Caption content |
|---|---|---|
| Crystal photograph before mounting | `docs/companion/assets/yzgo_crystal_before_mounting.jpg` | Crystal appearance, growth batch, scale bar if available. |
| Mounted sample photograph | `docs/companion/assets/yzgo_mounted_sample.jpg` | Mounting geometry, vertical-field orientation, neutron sample environment. |
| Orientation sketch | `docs/companion/assets/yzgo_orientation_sketch.png` | Crystallographic axes, scattering plane, field direction. |

Caption checklist:
- identify whether the image is the neutron sample or magnetization sample;
- include date/source if known;
- state the crystallographic orientation and field direction if visible;
- avoid over-interpreting color/shape features unless they are experimentally relevant.

## Run table source

The run inventory in this document was generated from the uploaded CSV:

`2026-06-16T14_36_51.277-04_00-SNS-CNCS-IPTS-37274-runs.csv`

Parsed top-level information:

| Field | Value |
|---|---|
| Sample name | `YbZn2GaO5` |
| Sample ID | `108367` |
| Run range | 663980–668461 |
| Number of rows/runs | 4482 |
| Total proton charge in CSV | 401.46 C |
| Start time range | 2025-10-17 18:11:27 to 2025-10-27 08:09:13 |
| Sample environment | `14T Vertical Field Magnet;Oxford Dil Fridge insert` |

## Core neutron datasets for the current analysis

The current mature analysis emphasizes the Ei = 4.65 meV and Ei = 3.32 meV CNCS datasets. Ei = 4.65 meV is the current co-fit neutron energy, while Ei = 3.32 meV is used as an important cross-check and background-comparison dataset because it has different kinematic coverage.

|   Ei (meV) |   T (K) |   B (T) |   runs |   first |   last |   charge (C) |   median C/run |   omega min |   omega max |   N omega |
|-----------:|--------:|--------:|-------:|--------:|-------:|-------------:|---------------:|------------:|------------:|----------:|
|       3.32 |    0.07 |       0 |    490 |  664977 | 665466 |        49.31 |           0.1  |        -190 |        54.5 |       490 |
|       3.32 |    0.07 |       9 |    491 |  665467 | 666102 |        49.14 |           0.1  |        -190 |        54.5 |       490 |
|       3.32 |    0.07 |      14 |    490 |  665957 | 666447 |        49.03 |           0.1  |        -190 |        54.5 |       490 |
|       3.32 |   20    |       0 |    490 |  664476 | 664965 |        49.01 |           0.1  |        -190 |        54.5 |       490 |
|       4.65 |    0.07 |       0 |    491 |  667971 | 668461 |        30.31 |           0.03 |        -190 |        54.5 |       490 |
|       4.65 |    0.07 |       9 |    490 |  667481 | 667970 |        49    |           0.1  |        -190 |        54.5 |       490 |
|       4.65 |    0.07 |      14 |    543 |  666448 | 667480 |        54.12 |           0.1  |        -190 |        54.5 |       490 |

Notes:
- The 3.32 meV dataset includes 0 T, 9 T, and 14 T at 0.07 K, plus a 20 K / 0 T dataset useful for background comparisons.
- The 4.65 meV dataset includes 0 T, 9 T, and 14 T at 0.07 K and is the current primary neutron co-fit dataset.
- The apparent 490 unique omega positions correspond to the nominal rotation scan from -190.0 deg to 54.5 deg in 0.5 deg steps. Some conditions have repeated or extra runs at selected omega values.

## Full run-condition summary

The full run CSV also includes alignment/check runs and additional incident-energy scans not currently emphasized in the co-fit.

| Ei (meV)   | T (K)   | B (T)   |   runs |   first |   last |   charge (C) |   median C/run | omega min   | omega max   |   N omega |
|:-----------|:--------|:--------|-------:|--------:|-------:|-------------:|---------------:|:------------|:------------|----------:|
| 3.32       | 0.07    | 0       |    490 |  664977 | 665466 |        49.31 |           0.1  | -190        | 54.5        |       490 |
| 3.32       | 0.07    | 9       |    491 |  665467 | 666102 |        49.14 |           0.1  | -190        | 54.5        |       490 |
| 3.32       | 0.07    | 14      |    490 |  665957 | 666447 |        49.03 |           0.1  | -190        | 54.5        |       490 |
| 3.32       | 20      | 0       |    490 |  664476 | 664965 |        49.01 |           0.1  | -190        | 54.5        |       490 |
| 4.65       | 0.07    | 0       |    491 |  667971 | 668461 |        30.31 |           0.03 | -190        | 54.5        |       490 |
| 4.65       | 0.07    | 9       |    490 |  667481 | 667970 |        49    |           0.1  | -190        | 54.5        |       490 |
| 4.65       | 0.07    | 14      |    543 |  666448 | 667480 |        54.12 |           0.1  | -190        | 54.5        |       490 |
| 4.65       | 12      | 0       |      1 |  664966 | 664966 |         1    |           1    | 0           | 0           |         1 |
| 4.75       | 12      | 0       |      1 |  664967 | 664967 |         1    |           1    | 0           | 0           |         1 |
| 4.85       | 12      | 0       |      1 |  664968 | 664968 |         1    |           1    | 0           | 0           |         1 |
| 4.9        | 12      | 0       |      1 |  664969 | 664969 |         1    |           1    | 0           | 0           |         1 |
| 5          | 12      | 0       |      1 |  664970 | 664970 |         1    |           1    | 0           | 0           |         1 |
| 9          | 12      | 0       |      1 |  664971 | 664971 |         1    |           1    | 0           | 0           |         1 |
| 9.25       | 12      | 0       |      1 |  664972 | 664972 |         1    |           1    | 0           | 0           |         1 |
| 9.5        | 12      | 0       |      1 |  664973 | 664973 |         1    |           1    | 0           | 0           |         1 |
| 12         | 12      | 0       |      1 |  664974 | 664974 |         1    |           1    | 0           | 0           |         1 |
| 12         |         | 0       |    491 |  663980 | 664470 |        10.56 |           0.02 | -190        | 54.5        |       490 |
| 13         | 0.07    | 14      |    490 |  666970 | 667459 |        49.01 |           0.1  | -190        | 54.5        |       490 |
| 13         | 12      | 0       |      1 |  664975 | 664975 |         1.01 |           1.01 | 0           | 0           |         1 |
| 25         | 12      | 0       |      1 |  664976 | 664976 |         1    |           1    | 0           | 0           |         1 |
|            |         |         |      5 |  664471 | 664475 |         0.96 |           0.03 |             |             |         0 |

## Reduced data products

Document the reduced data products with paths relative to the repository root.

Suggested table:

| Data product | Repo path | Used for | Generated by |
|---|---|---|---|
| CNCS 1D cuts | `data/neutron/CNCS_1d_scans/` | 1D neutron model comparison and co-fit | Mantid/reduction workflow, then repo scripts |
| CNCS 2D cuts/maps | `data/neutron/CNCS_2d_scans/` | 2D data-versus-model figures | Mantid/reduction workflow, then `scripts/plot_2d_data_vs_model.jl` |
| Magnetization CSV | `data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv` | M(H) comparison and co-fit | digitized magnetization curve |
| CIF | `data/cif/` | lattice/crystal information | crystallographic source file |

## Recommended figures for this section

1. Photograph of the YbZn2GaO5 crystal/sample.
2. Photograph or schematic of the mounted sample in the vertical-field geometry.
3. Run-condition histogram: incident energy versus field/temperature.
4. Omega coverage histogram for the primary 3.32 meV and 4.65 meV datasets.
5. Overview figure showing which datasets enter the primary co-fit and which serve as background/validation checks.

## Caveats and open items

- Confirm whether the neutron and magnetization samples are the same crystal or separate crystals.
- Add sample mass, shape, and orientation details.
- Add the final scattering-plane convention and sign convention used by the 1D and 2D data products.
- Decide whether additional Ei = 13 meV / 14 T or higher-energy check data should be described as supporting data or left in the run inventory only.




## 2. Neutron reduction, histogramming, and background subtraction

## Purpose

This section documents how the neutron inelastic data are converted into the 1D and 2D observables compared against the analytical model. The goal is to separate three issues that are easy to conflate:

1. experimental reduction and binning;
2. empirical background subtraction;
3. model deposition onto the same Q/E bins as the experiment.

## Primary reduced products

Current repo paths:

| Product | Path | Role |
|---|---|---|
| 1D cuts | `data/neutron/CNCS_1d_scans/` | Primary 9 T / 14 T co-fit data and 3.32/4.65 meV comparison cuts. |
| 2D cuts/maps | `data/neutron/CNCS_2d_scans/` | Data-versus-model maps generated from fixed best-fit parameters. |
| Background controls | `configs/background_compare_controls.toml` | Controls for comparing 4.65 meV and 3.32 meV background handling. |
| Background comparison script | `scripts/compare_1d_4p65_3p32_backgrounds.jl` | Repo-native driver for 3.32/4.65 meV comparison. |

## 1D scan convention

The co-fit currently uses three representative reciprocal-space cuts:

| qtag | Nominal position | Common interpretation |
|---|---|---|
| `0_1_0` | (0, 1, 0) | zone-center / equivalent Γ-like position depending on indexing convention |
| `0p33_0p33_0` | (1/3, 1/3, 0) | K-like position |
| `0p5_0_0` | (1/2, 0, 0) | M-like position |

Add here:
- precise HKL bin limits for each qtag;
- energy-bin width;
- normalization convention;
- whether intensities are per formula unit, per Yb, arbitrary units, or reduced-count units.

## Background handling

The current background-comparison workflow is motivated by the different kinematic coverage of the Ei = 4.65 meV and Ei = 3.32 meV datasets. The Ei = 4.65 meV data include the cuts used in the co-fit, while the Ei = 3.32 meV data include a 20 K / 0 T comparison dataset but do not cover the same Γ-point region.

Recommended narrative:

The 4.65 meV data are treated as the primary co-fit dataset. The 3.32 meV data are used to check whether spectral features and background assumptions are robust to a different incident energy and different kinematic coverage. For the 3.32 meV comparison, the available configurations include 0 T, 9 T, 14 T at 0.07 K and 0 T at 20 K. A conservative first-pass background estimate can be constructed from the minimum or low envelope across configurations, while the primary 4.65 meV analysis uses the established tail/background-subtracted mode.

## Model histogramming versus point evaluation

Document explicitly whether the model is evaluated at:
- the center of each experimental Q/E bin;
- a grid of points inside each experimental bin and averaged;
- a Monte Carlo sampling of the experimental bin volume.

This distinction matters because broad disorder distributions, finite instrument resolution, and rapidly dispersing modes can all change the apparent intensity and linewidth after integration over the experimental bin volume.

## Recommended figures

1. Raw and background-subtracted Ei = 4.65 meV 1D cuts at 9 T and 14 T.
2. Ei = 3.32 meV comparison cuts for 0 T, 9 T, 14 T, and 20 K / 0 T.
3. Overlay of 3.32 meV and 4.65 meV cuts after consistent background treatment.
4. Schematic showing model deposition into experimental Q/E bins.
5. One figure demonstrating sensitivity to background choice.

## Text to finalize later

- State the final canonical `data_mode` used in the co-fit.
- State which background comparison is qualitative versus which is used in fitted objective functions.
- Give a compact explanation of why the 3.32 meV data lack the same Γ coverage.




## 3. Analytical field-polarized neutron model with disorder

## Purpose

This section describes the fast analytical model used to compare field-polarized neutron spectra at 9 T and 14 T with the YbZn2GaO5 CNCS data. The model is the primary fitting backend in the current repository; Sunny.jl calculations are treated as a slower validation/comparison layer.

## Physical picture

The current model represents YbZn2GaO5 as an effective spin-1/2 triangular-lattice magnet in a field-polarized state. The dominant dispersive component is described by finite exchange interactions and disorder. A second, non-dispersive component is included phenomenologically to represent an effective independent-spin or weakly coupled contribution.

The model is intentionally pragmatic: it is designed to capture the field-polarized excitation energies, linewidths, and relative intensities in a way that can be evaluated quickly enough for repeated fitting.

## Hamiltonian-level parameters

Canonical parameter names are stored in `configs/best_fit_parameters.toml`.

| parameter                        |     value |
|:---------------------------------|----------:|
| gzz                              | 3.78458   |
| J1_meV                           | 0.229982  |
| J2_meV                           | 0.0087223 |
| sigma_gzz                        | 0.385924  |
| sigma_J                          | 0.239687  |
| gzz2                             | 2.99856   |
| sigma_gzz2                       | 0.843936  |
| gperp_ratio                      | 2.96772   |
| chi_vv_muB_per_T                 | 0.0782676 |
| second_kernel_relative_intensity | 0.157982  |
| neutron_global_scale             | 0.0448083 |
| magnetization_global_scale       | 0.418897  |

The microscopic/disorder parameters are:

| Parameter | Meaning | Units/convention |
|---|---|---|
| `J1_meV` | nearest-neighbor exchange scale | meV |
| `J2_meV` | next-nearest-neighbor exchange scale | meV |
| `gzz` | longitudinal g factor for dispersive component | dimensionless |
| `sigma_gzz` | disorder width for dispersive longitudinal g factor | dimensionless or fractional, confirm implementation convention |
| `sigma_J` | shared exchange-disorder width for J1 and J2 | fractional, confirm implementation convention |
| `gzz2` | longitudinal g factor for non-dispersive component | dimensionless |
| `sigma_gzz2` | disorder width for non-dispersive component | dimensionless or fractional, confirm implementation convention |
| `gperp_ratio` | effective transverse neutron-intensity ratio for flat/dispersive components | not a Hamiltonian g factor |
| `second_kernel_relative_intensity` | relative weight of non-dispersive component | dimensionless |
| `neutron_global_scale` | neutron intensity scale | arbitrary/current reduced-data units |

## Dispersive component

A compact form of the field-polarized model can be written as

```math
\mathcal{H} =
\sum_{\langle ij\rangle} J_1^{ij} \mathbf{S}_i\cdot\mathbf{S}_j
+
\sum_{\langle\langle ij\rangle\rangle} J_2^{ij} \mathbf{S}_i\cdot\mathbf{S}_j
-
\mu_B B \sum_i g_{zz}^i S_i^z .
```

For each disorder realization, the field-polarized spin-wave energy is evaluated from the exchange Fourier transform and the local Zeeman term. The implementation should be cited directly to the relevant functions in `src/` or the legacy co-fit script once this section is finalized.

Recommended final derivation steps to include:

1. define the triangular-lattice reciprocal-space convention;
2. define `J(q)` for J1 and J2 on that lattice;
3. give the polarized-mode energy expression used by the code;
4. state how the magnetic field enters through `gzz μB B`;
5. state how the neutron intensity factor is applied;
6. state the energy-resolution kernel used for comparison to experiment.

## Disorder averaging

The analytical model averages over disorder realizations. The current public configuration exposes a shared `sigma_J` for J1 and J2 and a `sigma_gzz` for the dispersive component. The intended conceptual model is:

```math
g_{zz}^i \sim P(g_{zz}; \bar{g}_{zz}, \sigma_{g}), \quad
J_n^{ij} \sim P(J_n; \bar{J}_n, \sigma_J).
```

Finalize by stating:
- distribution type: Gaussian, log-normal, truncated Gaussian, or other;
- whether `sigma_J` is absolute or fractional;
- how negative exchanges are handled if sampled;
- number of Monte Carlo samples used in production plots.

## Non-dispersive component

The second component is modeled as an effective S = 1/2 contribution with no dispersive exchange. Its excitation energy is controlled primarily by the Zeeman scale and its width by `sigma_gzz2`. In neutron observables, its relative intensity is controlled by `second_kernel_relative_intensity` and the effective transverse factor `gperp_ratio`.

Use language like:

The non-dispersive component should not be over-interpreted as a unique microscopic defect model at this stage. It is a compact representation of spectral weight that does not follow the primary dispersive branch within the current field-polarized analysis window.

## Generation of neutron observables

For each experimental cut, the model should reproduce the same observable as the data:

```math
I_{model}(E; Q_{bin}, B)
=
s_n \left[
I_{disp}(E; Q_{bin}, B)
+
r_2 g_{\perp,ratio}^2 I_{flat}(E; B)
\right],
```

where `s_n` is `neutron_global_scale` and `r2` is `second_kernel_relative_intensity`.

Document:
- whether the model is integrated over the full Q-bin volume;
- how the energy-resolution function is deposited;
- how negative background-subtracted points are treated in the objective function;
- whether the same relative flat/dispersive fraction is shared with magnetization.

## Recommended figures

1. Analytical dispersion at 9 T and 14 T along the relevant Q path.
2. Data/model overlays for the three 1D cuts at 9 T.
3. Data/model overlays for the three 1D cuts at 14 T.
4. Component-separated plots: dispersive, non-dispersive, and total intensity.
5. Sensitivity plot showing effects of `sigma_gzz`, `sigma_J`, and `sigma_gzz2`.

## Caveats

- The model is a field-polarized approximation and should not be used unchanged for zero-field correlations.
- Disorder is compressed into a small number of effective distributions.
- The non-dispersive component is phenomenological unless tied to an independent structural or spectroscopic assignment.
- Absolute neutron intensity normalization remains a key systematic uncertainty.




## 4. Analytical magnetization model with disorder

## Purpose

This section documents the model used to compare magnetization versus magnetic field with the same disorder and component structure used in the neutron analysis.

## Observable

The experimental observable is magnetization as a function of applied magnetic field, currently represented in the repository by

`data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv`

Add here:
- original measurement source;
- temperature;
- field direction;
- units before and after conversion;
- conversion to μB per Yb or per formula unit;
- whether any scale factor is fitted.

## Two-component structure

The magnetization model mirrors the neutron model conceptually:

```math
M_{total}(B,T) =
A_M \left[
(1-f_2) M_{disp}(B,T)
+
f_2 M_{flat}(B,T)
\right]
+
\chi_{vv} B .
```

Here:
- `M_disp` is the dispersive/mean-field component with exchange and g-factor disorder;
- `M_flat` is the effective non-dispersive S = 1/2 component;
- `f2` or the equivalent relative weight is tied to the same component structure used in the neutron model when the shared-fraction model is used;
- `chi_vv_muB_per_T` is a linear Van Vleck susceptibility term;
- `A_M` is `magnetization_global_scale`.

## Mean-field/disorder model

The mean-field magnetization should be described in terms of:
- local effective Zeeman field;
- exchange/mean-field correction;
- disorder averaging over `gzz` and exchange parameters;
- temperature-dependent spin-1/2 polarization;
- component averaging between dispersive and non-dispersive contributions.

A useful reference expression for an independent S = 1/2 moment is

```math
M(B,T;g) = \frac{g\mu_B}{2}\tanh\left(\frac{g\mu_B B}{2 k_B T}\right).
```

The final document should then state exactly how the dispersive mean-field part modifies this expression.

## Shared versus independent fractions

This repository has explored both separate and shared non-dispersive fractions. The current canonical co-fit naming indicates a shared-fraction workflow:

`results/fits/cofit_9T14T_shared_fraction`

The companion document should state clearly whether the neutron and magnetization components share the same relative fraction in the final reported parameter set.

## Recommended figures

1. Magnetization versus field with best-fit total model.
2. Component-separated magnetization: dispersive, non-dispersive, Van Vleck, total.
3. Overlay of low-temperature M(H) and an illustrative 2.5 K curve.
4. Sensitivity to `gzz2` and `sigma_gzz2` for the non-dispersive component.
5. Comparison of shared-fraction and independent-fraction models, if both are retained.

## Caveats

- If neutron and magnetization measurements used different crystals, component fractions need not be identical physically.
- The fitted `magnetization_global_scale` should be separated from microscopic moment physics in the narrative.
- Low-field behavior may be outside the field-polarized approximation, depending on the final fit window.




## 5. Sunny.jl validation layer

## Purpose

Sunny.jl calculations are included as an independent numerical validation layer for the analytical field-polarized model. They are not yet the primary fitting backend.

The companion document should keep the Sunny discussion separate from the analytical model so that readers understand which results are fitted and which are validation/comparison calculations.

## Current Sunny scripts

| Script | Purpose |
|---|---|
| `scripts/sunny_validate_magnetization_meanfield.jl` | Mean-field bridge/check against the analytical magnetization convention. |
| `scripts/sunny_validate_magnetization_largecell.jl` | Large-cell/disorder-style magnetization validation. |
| `scripts/sunny_validate_spinwave_kpm.jl` | 1D KPM spin-wave comparison to experimental cuts. |
| `scripts/sunny_plot_kpm_2d.jl` | 2D Sunny/KPM map generation. |
| `scripts/compare_analytical_sunny_dispersion.jl` | Analytical versus Sunny dispersion comparison. |
| `scripts/compare_analytical_vs_sunny_outputs.jl` | Higher-level comparison of analytical and Sunny outputs. |
| `scripts/check_sunny_form_factor.jl` | Checks built-in and fallback Yb3+ form-factor conventions. |
| `scripts/check_sunny_parameter_mapping.jl` | Checks mapping from canonical TOML parameters to Sunny conventions. |

## Current approximation choices

Document these explicitly:

1. The analytical model remains the fast fitting backend.
2. Sunny is used to validate the field-polarized results and to prepare for calculations outside the analytical approximation.
3. The flat component is represented as a zero-exchange or independent-spin Sunny component.
4. The KPM energy kernel is not identical to an experimental resolution convolution unless explicitly calibrated.
5. The Sunny transverse moment-tensor convention is a calculator/intensity gauge, not a fitted material `gperp`.
6. The physical `gperp_ratio` in the co-fit is an effective relative transverse-intensity factor outside the Hamiltonian.
7. The Yb3+ magnetic form factor should be handled consistently in 1D and 2D Sunny routes.

## Sunny intensity convention

For the current comparison, the total Sunny neutron model can be described schematically as

```math
I_{Sunny,total} =
s_{Sunny}
\left[
I_{disp,Sunny}
+
r_2 g_{\perp,ratio}^2 I_{flat,Sunny}
\right].
```

This mirrors the analytical component structure while allowing the Sunny comparison to have an independent intensity scale if desired.

## Recommended figures

1. Analytical dispersion versus Sunny dispersion at best-fit parameters.
2. Sunny KPM 1D comparison for the same 9 T / 14 T cuts used in the co-fit.
3. Sunny 2D maps at 9 T and 14 T.
4. Mean-field and large-cell magnetization comparison.
5. Form-factor diagnostic figure or table for Γ/K/M-like points.

## Open Sunny issues to track

- Decide whether the Sunny 1D workflow should histogram the full experimental Q-bin volume rather than evaluating only a representative path/grid.
- Decide how to map the KPM energy kernel onto the analytical resolution function.
- Record the Sunny.jl version used for every published figure.
- Record supercell size, random seed, number of disorder realizations, and KPM settings.




## 6. Reproducibility and provenance

## Purpose

This section maps each scientific claim or figure to the repository files that produce it. The companion document should make the analysis auditable: a reader should be able to find the input data, config file, script, and output file for each major result.

## Quick-start commands

From the repository root:

```bash
julia --project=. scripts/check_repo.jl
```

Primary co-fit:

```bash
julia --project=. scripts/run_cofit_9T14T.jl
```

Background comparison:

```bash
julia --project=. scripts/compare_1d_4p65_3p32_backgrounds.jl
```

2D data-versus-model plot:

```bash
julia --project=. scripts/plot_2d_data_vs_model.jl
```

Sunny validation examples:

```bash
julia --project=. scripts/sunny_validate_magnetization_meanfield.jl
julia --project=. scripts/sunny_validate_magnetization_largecell.jl
julia --project=. scripts/sunny_validate_spinwave_kpm.jl
```

## Core configs

| Config | Purpose |
|---|---|
| `configs/best_fit_parameters.toml` | Canonical parameter set / latest best-fit values. |
| `configs/cofit_controls.toml` | Observable selection, weights, sampling, and output name for the co-fit. |
| `configs/background_compare_controls.toml` | Controls for comparing Ei = 3.32 meV and Ei = 4.65 meV background handling. |
| `configs/feature_extraction_controls.toml` | Feature-extraction controls. |
| `configs/plot_2d_controls.toml` | Controls for 2D data-versus-model plotting. |
| `configs/sunny_validation_controls.toml` | Controls for Sunny validation calculations. |

## Draft figure provenance table

| Figure | Scientific purpose | Script | Config | Input data | Output path | Status |
|---|---|---|---|---|---|---|
| Crystal photograph | sample documentation | manual | n/a | photograph | `docs/companion/assets/...` | TODO |
| Run-condition histogram | experiment overview | TODO script | run CSV | run CSV | `results/figures/run_inventory/...` | TODO |
| 3.32/4.65 meV comparison | background validation | `scripts/compare_1d_4p65_3p32_backgrounds.jl` | `configs/background_compare_controls.toml` | `data/neutron/CNCS_1d_scans/` | `results/figures/...` | fill output |
| 9 T neutron co-fit overlays | primary fit | `scripts/run_cofit_9T14T.jl` | `configs/cofit_controls.toml`, `configs/best_fit_parameters.toml` | 1D cuts + magnetization CSV | `results/fits/cofit_9T14T_shared_fraction/` | fill output |
| 14 T neutron co-fit overlays | primary fit | `scripts/run_cofit_9T14T.jl` | same | same | same | fill output |
| Magnetization co-fit | primary fit | `scripts/run_cofit_9T14T.jl` | same | magnetization CSV | same | fill output |
| 2D data-versus-model | fixed model comparison | `scripts/plot_2d_data_vs_model.jl` | `configs/plot_2d_controls.toml`, `configs/best_fit_parameters.toml` | `data/neutron/CNCS_2d_scans/` | `results/figures/2d_data_vs_model/` | fill output |
| Sunny KPM 1D | validation | `scripts/sunny_validate_spinwave_kpm.jl` | `configs/sunny_validation_controls.toml` | 1D cuts | `results/figures/sunny_validation/` | active |
| Sunny magnetization | validation | `scripts/sunny_validate_magnetization_meanfield.jl` and/or `scripts/sunny_validate_magnetization_largecell.jl` | `configs/sunny_validation_controls.toml` | best-fit parameters | `results/figures/sunny_validation/` | active |

## Release checklist

Before tagging a release or pointing collaborators to the repository:

- [ ] Confirm `Project.toml` and `Manifest.toml` reproduce the analysis on a clean checkout.
- [ ] Run `scripts/check_repo.jl`.
- [ ] Run the primary co-fit script with `run_optimization = false` to regenerate figures from the canonical parameter set.
- [ ] Confirm all figure paths in the companion document exist.
- [ ] Add crystal/sample photographs and captions.
- [ ] Add final data-reduction notes for the CNCS cuts.
- [ ] Add a citation/acknowledgement section.
- [ ] Add a license and data-use statement if the repository is public.
- [ ] Create a GitHub release or archive with a DOI if the document is intended to be cited.

## Suggested citation block

Add final citation text here once the repository is ready.

Suggested placeholder:

> This repository contains analysis code and derived data products for the field-polarized neutron scattering and magnetization analysis of YbZn2GaO5. The canonical parameter set is stored in `configs/best_fit_parameters.toml`, and the companion documentation in `docs/companion/` describes the data, model assumptions, and reproducibility workflow.

