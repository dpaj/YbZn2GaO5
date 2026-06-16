# YbZn2GaO5 companion documentation starter

This directory is a starter companion-document layer for the `dpaj/YbZn2GaO5` GitHub repository. It is written as modular Markdown so that each section can mature independently while still forming a coherent reproducibility guide.

## Recommended structure

- `01_experiment_and_run_inventory.md` — sample photographs, neutron experiment inventory, run-condition summary, and data provenance.
- `02_neutron_reduction_backgrounds.md` — inelastic neutron cuts, histogramming, background subtraction, and 3.32/4.65 meV comparison.
- `03_analytical_field_polarized_neutron_model.md` — analytical field-polarized neutron model, disorder averaging, and experimental observable generation.
- `04_magnetization_model.md` — mean-field/disorder magnetization model and co-fit implementation notes.
- `05_sunny_validation.md` — Sunny.jl validation layer, current approximations, and comparison to the analytical model.
- `06_reproducibility_and_provenance.md` — script/config/figure mapping, run instructions, and release checklist.

## Current maturity map

| Topic | Status | Recommended next action |
|---|---|---|
| Neutron data collection and run histogram inventory | mature enough to document | Add final histogram figures and link raw/reduced data products. |
| Background subtraction and 3.32/4.65 meV comparison | mature enough to document | Add final comparison figure filenames and state which background mode is canonical. |
| Analytical field-polarized neutron model at 9 T and 14 T | mature enough to document | Add derivation equations checked against implementation. |
| Analytical model to neutron observables | mature enough to document | Document bin-center versus full histogram treatment explicitly. |
| Magnetization model and implementation | mature enough to document | Add units-conversion details and final M(H) component plots. |
| Sunny validation | active / partially mature | Keep separate from the primary analysis until approximation choices stabilize. |

## Tables generated from the attached run CSV

- `tables/neutron_run_summary_core_analysis.csv`
- `tables/neutron_run_summary_by_condition.csv`
- `tables/neutron_run_inventory_parsed.csv`

The parsed run CSV covers 4482 runs for sample `YbZn2GaO5`, run numbers 663980–668461, with total listed proton charge 401.46 C, from 2025-10-17 18:11:27 to 2025-10-27 08:09:13. The sample environment listed in the run table is `14T Vertical Field Magnet;Oxford Dil Fridge insert`.
