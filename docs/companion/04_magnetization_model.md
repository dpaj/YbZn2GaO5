# 4. Analytical magnetization model with disorder

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
