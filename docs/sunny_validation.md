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
