# YbZn2GaO5 neutron and magnetization co-fit

This repository contains analysis code and data products for modeling the field-polarized neutron scattering and magnetization of YbZn2GaO5.

The current fitting workflow uses a fast analytical polarized-phase model for the neutron and magnetization co-fit. A Sunny.jl version of the final simulation will be added as an independent validation/comparison calculation.

## Repository layout

- `configs/` — human-readable parameter and run-control files.
- `data/` — experimental inputs and crystallographic files.
- `src/` — reusable Julia functions.
- `scripts/` — runnable analysis and plotting scripts.
- `results/` — generated fit outputs, figures, and feature tables.
- `docs/` — notes on modeling choices, background subtraction, resolution, and reproducibility.