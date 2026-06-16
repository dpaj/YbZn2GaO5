# 2. Neutron reduction, histogramming, and background subtraction

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
