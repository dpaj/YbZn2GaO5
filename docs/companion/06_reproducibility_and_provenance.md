# 6. Reproducibility and provenance

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
