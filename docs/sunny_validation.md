# Sunny.jl validation workflow

This repo includes preliminary Sunny.jl validation scripts for the YbZn2GaO5 analytical co-fit model.

The fast analytical model remains the fitting backend.  Sunny.jl is used as a slower independent calculator to check the field-polarized results and, eventually, to compute zero-field magnetic correlations that are outside the simple analytical polarized-phase treatment.

## Current validation targets

The first Sunny workflows compare against the same canonical parameter set used by the analytical co-fit:

```text
configs/best_fit_parameters.toml
```

The validation scripts are:

```powershell
julia --project=. scripts/sunny_validate_magnetization_meanfield.jl
julia --project=. scripts/sunny_validate_magnetization_largecell.jl
julia --project=. scripts/sunny_validate_spinwave_kpm.jl
```

## Model convention

The canonical parameter set uses:

- `sigma_J`: one shared fractional exchange-disorder width for J1 and J2.
- `gperp_ratio`: effective flat/dispersive transverse neutron intensity ratio.
- `second_kernel_relative_intensity`: shared relative weight of the nondispersive / flat component.

The Sunny validation scripts should not expose legacy names such as `sigma_J1`, `sigma_J2`, `gperp`, or `gperp2` as fit parameters.  Any conversion to older helper-function conventions should occur only at an internal boundary.

## Total model

For the Sunny validation, the intended total model is:

```text
total = Sunny dispersive component + Sunny flat component
```

The flat component is represented as an independent-spin / zero-exchange Sunny component using `gzz2` and `sigma_gzz2`.  This is deliberately distinct from treating `gperp_ratio` as a physical g-factor.

## Status

These scripts are preliminary scaffolds.  The mean-field bridge is expected to run first.  The large-cell and KPM scripts may require Sunny API tuning depending on the installed Sunny version and the final choice of effective triangular-lattice builder.


## Sunny KPM 1D comparison layer

The KPM validation now loads the same Ei = 4.65 meV, T = 0.07 K, 9 T / 14 T 1D cuts targeted by the analytical co-fit.  The script histograms/interpolates the Sunny spectra onto the experimental energy grid before plotting.

Current preliminary convention:

```text
I_total(E) = neutron_global_scale * [I_disp_Sunny(E) + r2*gperp_ratio^2*I_flat_Sunny(E)]
```

The flat component is still a Sunny zero-exchange component; `gperp_ratio` is used only as the same effective relative neutron-intensity factor used in the analytical comparison, not as a physical Hamiltonian g-factor.

The helper script

```powershell
julia --project=. scripts/check_neutron_cut_loading.jl
```

prints the neutron cuts that will be used by the KPM validation and writes an inventory CSV.

## Sunny neutron intensity scale

The Sunny KPM neutron intensity is allowed to have an independent overall comparison scale from the analytical co-fit neutron scale.  The control

```toml
[kpm]
neutron_scale_mode = "least_squares"
neutron_scale_scope = "global"
```

fits one nonnegative global scale factor across the selected 1D cuts.  For quick visual diagnostics, set

```toml
neutron_scale_mode = "max_match"
```

to match the maximum Sunny model intensity to the maximum experimental intensity over the configured `neutron_scale_fit_window_meV`.

The analytical `neutron_global_scale` from `configs/best_fit_parameters.toml` is still available using

```toml
neutron_scale_mode = "best_fit"
```

and `manual_neutron_global_scale` is available with `neutron_scale_mode = "manual"`.

## Sunny transverse g convention for KPM neutron intensity

For the field-polarized H || c validation, the Sunny systems intentionally keep nonzero transverse moment-tensor components:

```toml
[common]
sunny_transverse_gxy = 1.0
```

This is a Sunny neutron-intensity gauge, not a fitted physical `gperp`.  The reason is that `ssf_perp(sys)` measures magnetic moment fluctuations; if the Sunny moment tensor has `gxx = gyy = 0`, the transverse inelastic neutron intensity can vanish numerically even when the longitudinal `gzz` and magnon energy are sensible.

The canonical physical comparison still keeps `gperp_ratio` outside the Hamiltonian as an effective flat/dispersive transverse intensity weight:

```text
I_total = neutron_scale * (I_disp_Sunny + r2 * gperp_ratio^2 * I_flat_Sunny)
```

Changing `sunny_transverse_gxy` should therefore be treated as a calculator/debug convention.  It should not be reported as a refined material g-factor.

## Sunny finite-size controls

Sunny system size is treated as an extrinsic calculator setting, not as a fitted physical parameter.  It is controlled in `configs/sunny_validation_controls.toml`.

For the large-cell magnetization validation:

```toml
[largecell]
system_size = [4, 4, 1]
dims = [4, 4, 1]
repeat_factor = [1, 1, 1]
```

For the KPM neutron validation:

```toml
[kpm]
system_size = [12, 12, 1]
dims = [3, 3, 1]
repeat_factor = [4, 4, 1]
```

`system_size` is the final finite Sunny supercell used for the calculation.  `dims` is the seed Sunny system constructed before `repeat_periodically`.  When `system_size` is present, the validation code derives

```text
repeat_factor = system_size ./ dims
```

and checks consistency with any explicitly supplied `repeat_factor`.

## Flat/dispersive S=1/2 fraction

The best-fit relative flat-component fraction is stored in the canonical best-fit file:

```toml
[neutron_extrinsic]
second_kernel_relative_intensity = 0.1579822309
```

The historical name `second_kernel_relative_intensity` means the best-fit relative weight of the nondispersive flat S=1/2 component compared with the dispersive component.  In the Sunny validation controls this is selected by:

```toml
[weights]
use_second_kernel_relative_intensity_from_best_fit = true
manual_second_kernel_relative_intensity = 0.1579822309
```

The CSV outputs now include both names where relevant:

```text
second_kernel_weight
flat_to_dispersive_fraction
```

For the neutron KPM comparison, the model still uses the effective transverse-intensity convention:

```text
I_total = neutron_scale * (I_disp_Sunny + flat_weight * I_flat_Sunny)
flat_weight = flat_to_dispersive_fraction * gperp_ratio^2
```

This keeps the fitted population/weight fraction distinct from the effective neutron matrix-element ratio.


## Parameter-mapping check

To verify that Sunny is using the same canonical parameter set as the analytical
co-fit, run:

```powershell
julia --project=. scripts/check_sunny_parameter_mapping.jl
```

This prints and writes:

```text
results/feature_tables/sunny_validation/sunny_parameter_mapping.csv
```

The important convention is:

```text
dispersive Sunny component:
    gzz, J1_meV, J2_meV, sigma_gzz, sigma_J

flat Sunny component:
    gzz2, sigma_gzz2, J1 = J2 = 0

neutron flat weight:
    flat_weight = second_kernel_relative_intensity * gperp_ratio^2
```

The Sunny transverse `gxy` is an intensity gauge required by `ssf_perp`; it is
not the fitted physical `gperp`.

## Sunny KPM 2D path map

A model-only 2D KPM map can be generated with:

```powershell
julia --project=. scripts/sunny_plot_kpm_2d.jl
```

The q path and field are controlled by:

```toml
[kpm_2d]
field_T = 9.0
qtags = ["0_1_0", "0p33_0p33_0", "0p5_0_0"]
n_per_segment = 41
neutron_scale_mode = "best_fit"
z_mode = "linear"
```

This writes:

```text
results/figures/sunny_validation/sunny_kpm_2d_path_9T.png
results/feature_tables/sunny_validation/sunny_kpm_2d_path_9T.csv
```

The 2D map is a Sunny model calculation only.  It does not yet overlay the
experimental 2D CNCS data or apply a fitted experimental 2D scale.
