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

## Yb3+ magnetic form factor in Sunny KPM

The Sunny KPM neutron validation now prefers Sunny's built-in magnetic form-factor machinery.  Both the 1D KPM comparison and the 2D KPM maps construct their neutron measurement through the same shared helper, so the form-factor setting applies to both routes.  The form factor is passed into the neutron measurement object, rather than being applied afterward as a manual `|f(Q)|^2` multiplier:

```julia
formfactors = [1 => FormFactor("Yb3")]
measure = ssf_perp(sys; formfactors)
swt = SpinWaveTheoryKPM(sys; measure, ...)
```

The control is:

```toml
[neutron_form_factor]
enabled = true
source = "sunny_builtin"
ion = "Yb3"
candidate_ions = ["Yb3", "Yb3+"]
on_error = "error"
```

This matches Sunny's documented usage pattern, where the tutorial example constructs the neutron measurement with `formfactors = [1 => FormFactor("Co2")]` and then calls `ssf_perp(sys; formfactors)`.  Keeping the form factor inside `ssf_perp` is cleaner than manually multiplying a final intensity map, and it avoids accidentally double counting the magnetic form factor.

The previous analytical/manual Yb3+ form-factor implementation is retained only as a diagnostic fallback:

```toml
source = "manual_yb3"
manual_include_j2 = true
manual_j2_coefficient = 0.75
manual_apply_as = "intensity_squared"
```

For normal Sunny KPM validation, leave `source = "sunny_builtin"`.  The physical reciprocal metric below is used only by the manual fallback and the diagnostic script's printed `|Q|` values:

```toml
[neutron_form_factor.lattice]
a_A = 3.376
c_A = 21.96
gamma_deg = 120.0
```

The helper script

```powershell
julia --project=. scripts/check_sunny_form_factor.jl
```

verifies that the Sunny built-in form-factor label can be constructed and prints the manual fallback `|Q|`, `f(Q)`, and `|f(Q)|^2` values at Γ, K, M, Γ1, and K1 for a quick sanity check.


### 1D / 2D consistency note

The 1D and 2D KPM paths should not call `ssf_perp(sys)` directly.  They should call the shared Sunny-validation helper that chooses either:

```julia
ssf_perp(sys; formfactors=[1 => FormFactor("Yb3")])
```

or, for diagnostics only, the no-form-factor or manual-fallback route.  The 1D CSV output now records the q center, `|Q|`, form-factor source, and form-factor diagnostic columns so that accidental divergence between the 1D and 2D routes is easier to catch.
