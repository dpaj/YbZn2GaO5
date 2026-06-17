# YbZn2GaO5 neutron and magnetization analysis

This repository contains analysis code, configuration files, reduced data products, documentation, and generated figures for modeling the field-polarized neutron scattering and magnetization of **YbZn2GaO5**.

The current workflow focuses on a fast analytical model for the field-polarized state at high magnetic field, with disorder in the effective magnetic parameters. The same modeling framework is used to compare against inelastic neutron scattering spectra and magnetization as a function of magnetic field. Sunny.jl-based calculations are maintained as an independent numerical validation/comparison layer.

## Current project status

The repository now contains several mostly mature pieces of the analysis:

- experimental data collection notes and neutron run inventory,
- reduced one-dimensional neutron inelastic spectra,
- background-subtraction routines,
- comparison of Ei = 4.65 meV and Ei = 3.32 meV neutron datasets,
- analytical field-polarized neutron model for 9 T and 14 T data,
- implementation of the analytical model to generate neutron observables,
- analytical mean-field/disorder model for magnetization,
- co-fit infrastructure for neutron and magnetization observables,
- Sunny.jl validation scripts and notes,
- companion documentation describing assumptions, approximations, and figure provenance.

The main analysis path currently uses the analytical polarized-state model because it is fast enough for fitting and parameter exploration. Sunny.jl calculations are used to check and interpret the analytical approximations, especially the consequences of disorder, finite-size effects, and numerical broadening.

## Repository layout

```text
configs/      Human-readable TOML parameter and run-control files
data/         Experimental inputs, reduced scans, magnetization data, and crystallographic files
docs/         Modeling notes, background-subtraction notes, Sunny validation notes, and companion docs
results/      Generated figures, fit outputs, feature tables, and validation products
scripts/      Runnable Julia scripts for analysis, plotting, fitting, and validation
src/          Reusable Julia source code for loading data, modeling, fitting, and plotting
```

The most important documentation layer is:

```text
docs/companion/
```

This folder contains modular companion documents covering the experiment, neutron reduction, analytical neutron model, magnetization model, Sunny validation, and reproducibility/provenance.

## Companion documentation

The companion documentation is intended to explain what the repository does and how the pieces connect scientifically.

Start here:

```text
docs/companion/README.md
```

Then see:

```text
docs/companion/01_experiment_and_run_inventory.md
docs/companion/02_neutron_reduction_backgrounds.md
docs/companion/03_analytical_field_polarized_neutron_model.md
docs/companion/04_magnetization_model.md
docs/companion/05_sunny_validation.md
docs/companion/06_reproducibility_and_provenance.md
```

A single combined draft is also available:

```text
docs/companion/companion_combined.md
```

The companion docs are still being refined, but they are intended to serve as the reproducibility and interpretation guide for the repository.

## Quick start

From the repository root, instantiate the Julia environment and check that the repository loads:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/check_repo.jl
```

To compare the Ei = 3.32 meV and Ei = 4.65 meV neutron 1D scans and background handling:

```powershell
julia --project=. scripts/compare_1d_4p65_3p32_backgrounds.jl
```

To run the main neutron/magnetization co-fit workflow:

```powershell
julia --project=. scripts/run_cofit_9T14T.jl
```

A faster smoke-test version is available as:

```powershell
julia --project=. scripts/run_cofit_9T14T_smoke.jl
```

To generate fixed-parameter comparison plots from the current best-fit model, see the plotting scripts in:

```text
scripts/
```

and the run-control files in:

```text
configs/
```

## Data included

The repository includes reduced data products used by the analysis, including:

- one-dimensional neutron scattering cuts,
- magnetization versus magnetic field data,
- crystallographic input files,
- configuration files defining fitting and plotting controls,
- derived feature tables and generated figures.

Large raw neutron event files are not expected to be stored directly in this repository. The companion documentation records the run inventory and data provenance needed to connect the reduced data products to the original experiment.

## Modeling overview

### Neutron scattering model

The high-field inelastic neutron data are modeled using an analytical field-polarized spin-wave-like treatment for an effective spin-1/2 triangular-lattice model. The model includes:

- nearest-neighbor and further-neighbor exchange parameters,
- anisotropic g-factor parameters,
- disorder in gzz,
- disorder in exchange,
- a dispersive magnetic excitation component,
- an additional non-dispersive or weakly dispersive component,
- broadening and binning procedures to compare calculated spectra with experimental cuts.

The analytical model is used because it is computationally efficient enough for fitting and systematic parameter exploration.

### Magnetization model

The magnetization is modeled using a mean-field/disorder treatment with the same physical parameter set where possible. The magnetization model includes:

- disorder-averaged field response,
- Van Vleck susceptibility,
- dispersive and non-dispersive components,
- finite-temperature magnetization curves,
- comparison to measured M(H) data.

### Sunny.jl validation

Sunny.jl scripts are included as an independent validation and comparison layer. These calculations are used to test the analytical approximations and to understand where numerical disorder, finite-size effects, kernel broadening, and histogramming choices may matter.

At present, the analytical model is the primary fitting engine, while Sunny.jl is used for validation, interpretation, and future extensions. Ultimately, if acceleration is sufficient testing of optimization with Sunny.jl as the engine.

## Reproducibility philosophy

The goal of this repository is not only to store scripts, but to make the analysis reproducible and auditable.

Important reproducibility information is organized around:

- input data files,
- configuration files,
- analysis scripts,
- generated figures,
- model parameters,
- figure provenance,
- known approximations and limitations.

The companion documentation should make it possible to answer:

1. What data were used?
2. What model was fit?
3. What approximations were made?
4. Which script generated each figure?
5. Which configuration file and parameter set were used?
6. What parts of the analysis are mature, and what parts are still under active development?

## Suggested reading order

For a scientific overview:

1. `docs/companion/README.md`
2. `docs/companion/01_experiment_and_run_inventory.md`
3. `docs/companion/02_neutron_reduction_backgrounds.md`
4. `docs/companion/03_analytical_field_polarized_neutron_model.md`
5. `docs/companion/04_magnetization_model.md`
6. `docs/companion/05_sunny_validation.md`

For reproducing figures or checking scripts:

1. `docs/companion/06_reproducibility_and_provenance.md`
2. `configs/`
3. `scripts/`
4. `results/`

## Known active-development areas

The following items are still evolving:

- final choice of canonical background subtraction for some neutron cuts,
- quantitative comparison of analytical histogramming versus Sunny/KPM-style calculations,
        Analytical sigma_J:
            distribution of clean dispersions from different effective regions

        Sunny bond sigma_J:
            one disordered Hamiltonian with many exchange values coupled together
            
- absolute neutron intensity scale and scale-factor interpretation (some weirdness for the phase fractions),
- finite-size and disorder-realization effects in Sunny.jl,
- GPU acceleration of LSWT KPM calcs,
- AI utilization and more global optimization,
- Optimization of parameters based on Sunny,
- Finite temperature for the Sunny.jl magnetization large cell model,
- utilization of Sunny.jl to extend high-field models to zero-field,
- final figure provenance table for all publication-quality outputs.

## Citation and acknowledgement

...

## License

...
