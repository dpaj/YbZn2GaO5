# NHMFL SCM1 AC susceptibility, July/August 2025

48 sweeps from the SCM1 dilution-refrigerator insert at the NHMFL, plus the run log
`SCM1_July2025.xlsx`. Raw instrument text, tab-separated, one file per sweep.

## READ THE RUN LOG BEFORE THE COLUMN NAMES

The probe carried **three pickup coils and only one held our sample.** The column names
(`B1x1`, `T1x1`, `T3x1`, ...) name coil positions, not samples, and the mapping is only in
the xlsx:

| coil | contents |
|---|---|
| `T3` | **YbZn2GaO5, "para" = B ∥ c — OURS, and the only one** |
| `T1` | LuCu(OH)Br — a *different compound*, another group sharing the probe |
| `B1` | **not listed in the log**, yet carries the largest signal (6x T3) |

So **there is no perpendicular AC measurement in this dataset.** An early analysis here
assumed `B1` was the perpendicular crystal. That was wrong, and self-refuting: `B1`'s
in-phase channel crosses zero and goes negative above ~12.4 T, which no sample
susceptibility can do, and its phase rotates 79 degrees over 0-12 T where `T3`'s moves 4.

## Lock-ins

Three SR830s at 991, 313 and 137 Hz, but **only 991 Hz drives** in the field sweeps, so the
`x2`/`x3`/`y2`/`y3` columns are noise at 1e-9 -- three orders below signal. Use `x1`,`y1`
only, and prefer the magnitude `sqrt(x^2+y^2)`, which is free of the lock-in phase
convention. Frequencies rotate between runs; the log records which.

## Temperature

`Tmc` is the mixing-chamber thermometer and is the one to use. **The run log's milliamp
figures are heater currents, not fields**, so runs that look like repeats are not:

| runs | sweep | Tmc |
|---|---|---|
| 015 / 016 | 0 to 18 T / 18 to 0 T | **20 mK** |
| 047 / 048 | 0 to 18 T / 18 to 0 T | **450-500 mK** (048 drifts to 0.60 K) |

## What this dataset can and cannot support

**Cannot: chi' itself.** `T3`'s magnitude is flat from ~1 to ~10 T and falls to ~40% at 18 T,
while the DC dM/dH falls to ~5% of its 1 T value by 10 T -- a factor ~8 disagreement, in the
wrong direction for temperature to explain, since cooling sharpens a saturation rather than
flattening one. Two background models were tried and both fail: a constant complex coil
offset leaves the flat region intact (6.77e-7 at 2 T vs 6.63e-7 at 8 T), and differencing
the 20 mK and 450 mK runs gives a residual that *grows* to 7.8e-8 V at 14-18 T, largest
exactly where the sample is most saturated.

**Can: bound sharp features.** A smooth instrumental background cannot create or cancel a
sharp feature, so the absence of any dip, step, spike or up/down hysteresis anywhere in
0-18 T does bound first-order transitions and magnetisation-plateau edges. In both sweep
pairs the minimum of chi' above 2 T sits at the 18 T endpoint -- there is no interior dip.

**To make it quantitative, three things are needed from the experiment**, not from analysis:
an empty-coil run at matching field and temperature, the coil constants, and an account of
what coil `B1` held.

Analysis: `scripts/plot_ac_susceptibility_plateau_test.jl`.
