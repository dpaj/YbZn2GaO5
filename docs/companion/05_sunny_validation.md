# 5. Sunny.jl validation layer

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
