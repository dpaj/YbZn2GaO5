# Background handling for YbZn2GaO5 CNCS data

This note documents the background-handling checks used for the YbZn2GaO5 neutron analysis.

## Motivation

The Ei = 4.65 meV data and Ei = 3.32 meV data have different kinematic coverage. The Ei = 4.65 meV data include the relevant low-temperature field-dependent cuts used in the co-fit. The Ei = 3.32 meV data provide additional comparison data, including 20 K at 0 T, but do not contain the same gamma-point coverage because of neutron kinematic constraints.

## Current comparison workflow

The script

```powershell
julia --project=. scripts/compare_1d_4p65_3p32_backgrounds.jl