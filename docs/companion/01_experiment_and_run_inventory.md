# 1. Experiment and run inventory

## Purpose

This section documents the experimental data products that underlie the YbZn2GaO5 neutron and magnetization analysis. The goal is to make it clear which measurements were performed, which subset is used in the current modeling workflow, and how the run metadata maps onto reduced data files and analysis scripts.

## Sample and crystal photographs

Add crystal/sample photographs here.

Suggested files:

| Placeholder | Suggested path | Caption content |
|---|---|---|
| Crystal photograph before mounting | `docs/companion/assets/yzgo_crystal_before_mounting.jpg` | Crystal appearance, growth batch, scale bar if available. |
| Mounted sample photograph | `docs/companion/assets/yzgo_mounted_sample.jpg` | Mounting geometry, vertical-field orientation, neutron sample environment. |
| Orientation sketch | `docs/companion/assets/yzgo_orientation_sketch.png` | Crystallographic axes, scattering plane, field direction. |

Caption checklist:
- identify whether the image is the neutron sample or magnetization sample;
- include date/source if known;
- state the crystallographic orientation and field direction if visible;
- avoid over-interpreting color/shape features unless they are experimentally relevant.

## Run table source

The run inventory in this document was generated from the uploaded CSV:

`2026-06-16T14_36_51.277-04_00-SNS-CNCS-IPTS-37274-runs.csv`

Parsed top-level information:

| Field | Value |
|---|---|
| Sample name | `YbZn2GaO5` |
| Sample ID | `108367` |
| Run range | 663980–668461 |
| Number of rows/runs | 4482 |
| Total proton charge in CSV | 401.46 C |
| Start time range | 2025-10-17 18:11:27 to 2025-10-27 08:09:13 |
| Sample environment | `14T Vertical Field Magnet;Oxford Dil Fridge insert` |

## Core neutron datasets for the current analysis

The current mature analysis emphasizes the Ei = 4.65 meV and Ei = 3.32 meV CNCS datasets. Ei = 4.65 meV is the current co-fit neutron energy, while Ei = 3.32 meV is used as an important cross-check and background-comparison dataset because it has different kinematic coverage.

|   Ei (meV) |   T (K) |   B (T) |   runs |   first |   last |   charge (C) |   median C/run |   omega min |   omega max |   N omega |
|-----------:|--------:|--------:|-------:|--------:|-------:|-------------:|---------------:|------------:|------------:|----------:|
|       3.32 |    0.07 |       0 |    490 |  664977 | 665466 |        49.31 |           0.1  |        -190 |        54.5 |       490 |
|       3.32 |    0.07 |       9 |    491 |  665467 | 666102 |        49.14 |           0.1  |        -190 |        54.5 |       490 |
|       3.32 |    0.07 |      14 |    490 |  665957 | 666447 |        49.03 |           0.1  |        -190 |        54.5 |       490 |
|       3.32 |   20    |       0 |    490 |  664476 | 664965 |        49.01 |           0.1  |        -190 |        54.5 |       490 |
|       4.65 |    0.07 |       0 |    491 |  667971 | 668461 |        30.31 |           0.03 |        -190 |        54.5 |       490 |
|       4.65 |    0.07 |       9 |    490 |  667481 | 667970 |        49    |           0.1  |        -190 |        54.5 |       490 |
|       4.65 |    0.07 |      14 |    543 |  666448 | 667480 |        54.12 |           0.1  |        -190 |        54.5 |       490 |

Notes:
- The 3.32 meV dataset includes 0 T, 9 T, and 14 T at 0.07 K, plus a 20 K / 0 T dataset useful for background comparisons.
- The 4.65 meV dataset includes 0 T, 9 T, and 14 T at 0.07 K and is the current primary neutron co-fit dataset.
- The apparent 490 unique omega positions correspond to the nominal rotation scan from -190.0 deg to 54.5 deg in 0.5 deg steps. Some conditions have repeated or extra runs at selected omega values.

## Full run-condition summary

The full run CSV also includes alignment/check runs and additional incident-energy scans not currently emphasized in the co-fit.

| Ei (meV)   | T (K)   | B (T)   |   runs |   first |   last |   charge (C) |   median C/run | omega min   | omega max   |   N omega |
|:-----------|:--------|:--------|-------:|--------:|-------:|-------------:|---------------:|:------------|:------------|----------:|
| 3.32       | 0.07    | 0       |    490 |  664977 | 665466 |        49.31 |           0.1  | -190        | 54.5        |       490 |
| 3.32       | 0.07    | 9       |    491 |  665467 | 666102 |        49.14 |           0.1  | -190        | 54.5        |       490 |
| 3.32       | 0.07    | 14      |    490 |  665957 | 666447 |        49.03 |           0.1  | -190        | 54.5        |       490 |
| 3.32       | 20      | 0       |    490 |  664476 | 664965 |        49.01 |           0.1  | -190        | 54.5        |       490 |
| 4.65       | 0.07    | 0       |    491 |  667971 | 668461 |        30.31 |           0.03 | -190        | 54.5        |       490 |
| 4.65       | 0.07    | 9       |    490 |  667481 | 667970 |        49    |           0.1  | -190        | 54.5        |       490 |
| 4.65       | 0.07    | 14      |    543 |  666448 | 667480 |        54.12 |           0.1  | -190        | 54.5        |       490 |
| 4.65       | 12      | 0       |      1 |  664966 | 664966 |         1    |           1    | 0           | 0           |         1 |
| 4.75       | 12      | 0       |      1 |  664967 | 664967 |         1    |           1    | 0           | 0           |         1 |
| 4.85       | 12      | 0       |      1 |  664968 | 664968 |         1    |           1    | 0           | 0           |         1 |
| 4.9        | 12      | 0       |      1 |  664969 | 664969 |         1    |           1    | 0           | 0           |         1 |
| 5          | 12      | 0       |      1 |  664970 | 664970 |         1    |           1    | 0           | 0           |         1 |
| 9          | 12      | 0       |      1 |  664971 | 664971 |         1    |           1    | 0           | 0           |         1 |
| 9.25       | 12      | 0       |      1 |  664972 | 664972 |         1    |           1    | 0           | 0           |         1 |
| 9.5        | 12      | 0       |      1 |  664973 | 664973 |         1    |           1    | 0           | 0           |         1 |
| 12         | 12      | 0       |      1 |  664974 | 664974 |         1    |           1    | 0           | 0           |         1 |
| 12         |         | 0       |    491 |  663980 | 664470 |        10.56 |           0.02 | -190        | 54.5        |       490 |
| 13         | 0.07    | 14      |    490 |  666970 | 667459 |        49.01 |           0.1  | -190        | 54.5        |       490 |
| 13         | 12      | 0       |      1 |  664975 | 664975 |         1.01 |           1.01 | 0           | 0           |         1 |
| 25         | 12      | 0       |      1 |  664976 | 664976 |         1    |           1    | 0           | 0           |         1 |
|            |         |         |      5 |  664471 | 664475 |         0.96 |           0.03 |             |             |         0 |

## Reduced data products

Document the reduced data products with paths relative to the repository root.

Suggested table:

| Data product | Repo path | Used for | Generated by |
|---|---|---|---|
| CNCS 1D cuts | `data/neutron/CNCS_1d_scans/` | 1D neutron model comparison and co-fit | Mantid/reduction workflow, then repo scripts |
| CNCS 2D cuts/maps | `data/neutron/CNCS_2d_scans/` | 2D data-versus-model figures | Mantid/reduction workflow, then `scripts/plot_2d_data_vs_model.jl` |
| Magnetization CSV | `data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv` | M(H) comparison and co-fit | digitized magnetization curve |
| CIF | `data/cif/` | lattice/crystal information | crystallographic source file |

## Recommended figures for this section

1. Photograph of the YbZn2GaO5 crystal/sample.
2. Photograph or schematic of the mounted sample in the vertical-field geometry.
3. Run-condition histogram: incident energy versus field/temperature.
4. Omega coverage histogram for the primary 3.32 meV and 4.65 meV datasets.
5. Overview figure showing which datasets enter the primary co-fit and which serve as background/validation checks.

## Caveats and open items

- Confirm whether the neutron and magnetization samples are the same crystal or separate crystals.
- Add sample mass, shape, and orientation details.
- Add the final scattering-plane convention and sign convention used by the 1D and 2D data products.
- Decide whether additional Ei = 13 meV / 14 T or higher-energy check data should be described as supporting data or left in the run inventory only.
