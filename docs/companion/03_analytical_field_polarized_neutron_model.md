# 3. Analytical field-polarized neutron model with disorder

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
