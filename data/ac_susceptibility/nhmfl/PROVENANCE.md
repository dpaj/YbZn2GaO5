# NHMFL SCM1 AC susceptibility, July/August 2025

48 sweeps from the SCM1 dilution-refrigerator insert at the NHMFL, plus the run log
`SCM1_July2025.xlsx`. Raw instrument text, tab-separated, one file per sweep.

## The coil key — read it from the xlsx, do not infer it from column names

The column names (`B1x1`, `T1x1`, `T3x1`, ...) name **coil positions, not samples**. The probe
carried three coils and the mapping is in the spreadsheet's **columns E, F, G, rows 1-2**:

| coil | contents |
|---|---|
| `B1` | **YbZn2GaO5, "perp" = B ⟂ c** |
| `T1` | LuCu(OH)Br — a *different compound*, another group sharing the probe |
| `T3` | **YbZn2GaO5, "para" = B ∥ c** |

**Both orientations of our sample are present.** An earlier analysis in this repo claimed there
was no perpendicular measurement; that was a spreadsheet-parsing error on our side (a
self-closing empty cell swallowed the neighbouring one, so `B1`'s entry was read as a bare
integer and mistaken for a row label). It also argued that `B1` could not be a sample because
its in-phase channel goes negative above ~12.4 T. **That argument was wrong**: a raw lock-in
quadrature contains the coil's own mutual-inductance background, which is large here and can
dominate, so either quadrature may be negative. Nothing about a negative `x1` rules out a sample.

Below rows 1-2, columns E/F/G are reused for the **per-run lock-in drive settings** (voltage and
frequency), matching the `Set SR830 ( GPIB0::n )` lines in each data file's own sequence history.

## Lock-ins

Three SR830s at 991, 313 and 137 Hz (GPIB 1, 4, 7). In the field sweeps **only the 991 Hz drive
produces signal**: `x1`/`y1` are ~1e-6 to 1e-7 while `x2`/`x3`/`y2`/`y3` sit at 1e-9, three orders
down. Use `x1`,`y1`. Work in the **complex plane** (`x1 + i*y1`) rather than with either quadrature
alone — the coil background and the sample response have different phases, so a single quadrature
mixes them, and the magnitude alone can be non-monotonic even when the sample response is not.

## Temperature

`Tmc` is the mixing-chamber thermometer and is the one to use. **The run log's milliamp figures are
heater currents, not fields**, so runs that look like repeats are not:

| runs | sweep | Tmc |
|---|---|---|
| 015 / 016 | 0 to 18 T / 18 to 0 T | **20 mK** |
| 047 / 048 | 0 to 18 T / 18 to 0 T | **450-500 mK** (048 drifts to 0.60 K) |

## The instrumental background, and the in-situ reference that partly beats it

There is a large field-dependent instrumental term. Referencing each channel to its own 18 T value
does not remove it — all three coils, *including the different compound*, then show a common
near-linear ramp to zero, which no set of three different samples would produce.

**`T1` is the useful handle.** LuCu(OH)Br sits on the same probe in the same magnet at the same
temperature and does not saturate over this range, so its field dependence is a template for the
instrumental drift. Subtracting a fitted multiple of `T1`'s complex field dependence from each
YZGO channel gives, each normalised to its own 1 T value:

| B (T) | perp AC, 20 mK | perp DC dM/dH, 2.5 K | para AC, 20 mK | para DC dM/dH, 2.5 K |
|---|---|---|---|---|
| 1 | 1.000 | 1.000 | 1.000 | 1.000 |
| 2 | 0.913 | 0.922 | 0.844 | 0.882 |
| 3 | 0.817 | 0.820 | 0.772 | 0.752 |
| 5 | 0.345 | 0.691 | 0.724 | 0.452 |
| 8 | 0.109 | 0.371 | 0.590 | 0.125 |
| 12 | 0.056 | — | 0.277 | 0.034 |

**The perpendicular channel works.** It agrees with DC to ~1% at 1-3 T and then saturates *faster*,
which is the correct direction for 20 mK against 2.5 K — cooling sharpens a saturation, and the DC
comparison is at a temperature 125x higher. **The parallel channel does not**: it stays high where
the DC data say the sample is saturated, and that is still unexplained.

Two honest caveats. The scale on the `T1` template is **fitted on 10-18 T assuming the sample is
saturated there**, so the high-field end of the perp agreement is partly circular; the 1-8 T
comparison is not. And `T1` is a different sample in a different coil, so it is a template for the
drift's *shape*, not a calibrated background. (The fitted scales are also large and of opposite sign
for the two channels, which is a further reason to treat this as provisional.)

## What this dataset supports

- **A sharp-feature bound, robustly.** A smooth background cannot create or cancel a sharp feature,
  so the absence of any dip, step, spike or up/down hysteresis anywhere in 0-18 T bounds
  first-order transitions and magnetisation-plateau edges. In both sweep pairs the minimum above
  2 T sits at the 18 T endpoint — no interior dip. This conclusion needs no calibration at all.
- **A perpendicular chi'(H) shape at 20 mK**, usable at 1-8 T with the `T1` reference, consistent
  with DC. This is a genuinely independent observable at fields past the 14 T DC and neutron data.
- **Not** an absolute chi', and **not** the parallel orientation.

**Still worth asking the experiment for**: an empty-coil run at matching field and temperature, and
the coil constants. Those would replace the fitted `T1` template with a measured background and
would settle the parallel channel.

Analysis: `scripts/plot_ac_susceptibility_plateau_test.jl`.
