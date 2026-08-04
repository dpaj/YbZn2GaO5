# Primary M(H) at 0.4 K, B parallel to c

**Put `YbZnGaO_Bpara001_12.2mg_MvH_0.4K_10042025.dat` here, unrenamed.**

This is the **primary measurement** replacing
`data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv`, which was digitised from a
figure. Same temperature regime as everything the M(H) work has been fitted against, so this is a
direct upgrade rather than a new observable.

Sample mass **12.2 mg** (from the filename), **B || c**, nominal **0.4 K**.

## Why this matters more than the 2.5 K files

The digitised curve is demonstrably mis-normalised. At 6.975 T the 2.5 K PPMS measurement gives
**1.653 uB/Yb** for B || c, while the digitised 0.4 K curve reads **~1.12** at the same field --
and 0.4 K must give MORE magnetisation than 2.5 K, not 32% less. So the digitised data are low by
at least a factor 1.48. This file removes that problem entirely by giving absolute uB/Yb from a
known mass.

It should also remove the "scar" in the digitised trace, which was an artefact of reading points
off a printed figure.

## What it will change, and why the switch must be deliberate

Every M(H) number in this repo currently derives from the digitised curve: `A_M`, `chi_vv`, the
rms values, and the `B_sat ~ J1/gzz` determination. Switching the input changes all of them.

`A_M` was always profiled out, so a pure scale error never affected any SHAPE fit -- but it does
mean the fitted `A_M = 0.373` was never interpretable as physics.

The interesting one is **open thread 5**: the fitted linear term came out at
`chi_vv = 0.19 uB/T` against a crystal-field allowance of `0.0171 +- 0.0007`, eleven times too
large, and the M(H) fit leaned on it heavily (forcing it to zero degraded rms from 0.0107 to
0.0652). If that was partly an artefact of a distorted digitised shape rather than a real physical
term, this file will show it. If the excess survives against primary data, it is real and needs a
physical explanation.

## Metadata still needed here

- [ ] **Which instrument.** 0.4 K is below a standard PPMS base of 1.8 K, so this needed a He-3 or
      dilution insert. Record which, because it bears on thermometry accuracy at the bottom of the
      range -- and the M(H) protocol treats 0.42 K as effectively T = 0.
- [ ] **Whether this is the same crystal** as the neutron sample, and as the 4.81 mg B || c crystal
      used for the 2.5 K measurement. Three masses have appeared (12.2, 4.81, 7.7 mg) and it is not
      obvious how many distinct crystals that represents.
- [ ] **Background / holder subtraction status**, and demagnetisation correction status with crystal
      shape if uncorrected.
- [ ] **Field range and step**, and whether it is a single ramp or a loop.

## Handling plan

Conversion is `mu_B/Yb = emu / (5585 * moles)` with `moles = mass_g / 453.53`, using
`M(YbZn2GaO5) = 453.53 g/mol` and one Yb per formula unit.

Rather than teach the loader a new format, a tracked conversion script writes a derived CSV in the
same two-column shape `sv_read_magnetization_csv` already expects. That keeps the raw file, the
conversion, and the derived input all auditable and re-runnable, and makes the config switch a
one-line change.

**Do not point the configs at the new data until the digitised-versus-primary overplot has been
looked at.** The scale factor and any shape difference should be seen explicitly, not discovered
as a silent change in every M(H) result.
