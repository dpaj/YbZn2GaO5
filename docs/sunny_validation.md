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

## Magnetization normalization and scale convention

The Sunny magnetization validation now follows the analytical co-fit comparison
convention explicitly.

The two magnetic components are first combined using the shared second-kernel
relative weight

```text
Mmag = (Mdisp + r2*Mflat)/(1+r2)
```

when `normalize_second_kernel_weight = true`.

The Van Vleck-like term is then added to form the unscaled comparison curve:

```text
Mcombo_unscaled = Mmag + chi_vv * B
```

Finally, the analytical co-fit magnetization scale from

```text
configs/best_fit_parameters.toml

[magnetization_extrinsic]
magnetization_global_scale = ...
```

is applied to the full unscaled curve:

```text
Mtotal = magnetization_global_scale * Mcombo_unscaled
```

This matches the analytical co-fit CSV convention, where `M_disp_scaled`,
`M_vv_scaled`, and `M_nondispersive_scaled` are all scaled by the same
`magnetization_scale`.

The Sunny output CSVs retain raw, unscaled, and scaled component columns so that
spin-vs-moment sign/scale conventions remain diagnosable.

The large-cell `minimize_energy!` workflow also has an explicit `moment_sign`
option under `[largecell]` because Sunny spin and magnetic-moment sign
conventions can differ from the experimental plotting convention.
