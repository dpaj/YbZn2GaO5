# AC susceptibility vs field, dilution-fridge temperatures, NHMFL (Tallahassee)

**Put the raw scan files and the Excel run log here, unrenamed.** The spreadsheet can stay as
`.xlsx` -- `openpyxl` and `pandas` are both available, so no CSV export is needed.

The run log contains YZGO in both **B || c** and **B perpendicular to c**, and also an unrelated
**Lu-containing compound which is to be ignored**. Recording which is which here means nobody has
to re-derive that mapping from the spreadsheet later.

## What is actually measured, and why it matters

**The measured quantity is dM/dH, not M.** M is a *derived* quantity obtained by integrating, so:

- the field derivative is the primary, least-processed observable and should be compared directly
  against the model's dM/dH rather than against an integrated M;
- integration introduces an arbitrary constant and accumulates any low-field artefact across the
  whole range, so an integrated M inherits errors that the raw derivative does not;
- **there are no absolute units on the magnetisation.** Only shapes are comparable, with a free
  scale -- which is already how every neutron and M(H) comparison in this repo works, so it costs
  nothing.

## Why this dataset may be the sharpest test of the central claim

A magnetisation **plateau appears in dM/dH as a dip toward zero**, and a field-induced transition
as a step or a spike. The published model requires structured non-k=0 phases with a clear plateau,
and this dataset is field-swept dM/dH at **dilution temperature and high field** -- exactly the
regime where those features are predicted.

A plateau is a **scale-free shape feature**, so the absence of units does not weaken this use at
all, and the higher field range is a genuine advantage over the 14 T PPMS data. This complements
the high-field neutron diffraction (which sees only k = 0) on the same question from an independent
probe.

So the primary intended use is **testing for the absence of the predicted phases**, not refining
parameters. Whether it becomes a fit input is a later and separate decision -- the artefacts and
missing units argue for using it as a qualitative discriminator first.

## Metadata needed here

- [ ] **Run-to-file mapping**: which scan files are YZGO B || c, which are YZGO B perp c, and which
      belong to the Lu compound and should be skipped. (The spreadsheet has this; recording the
      conclusion here avoids re-deriving it.)
- [ ] **Temperature(s)** of each scan, and whether thermometry is reliable under a field sweep at
      dilution temperatures.
- [ ] **Field range and sweep rate**, and whether up and down sweeps are both present -- sweep-rate
      dependence and hysteresis are worth knowing about separately from physics.
- [ ] **Whether the files contain raw dM/dH, an already-integrated M, or both.**
- [ ] **Known artefacts** to expect and their field locations, since the point of the dataset is
      finding shape features and an artefact would masquerade as one. This is the same hazard the
      Ei = 4.65 magnet background created for the neutron cuts.
- [ ] **Sample identity / mass** if known, and whether these are the same crystals as the PPMS and
      neutron measurements. All crystals so far come from the same Haidong Zhou group growth (Aya
      Rutherford), the same growth as the neutron crystal, so cross-comparison is legitimate.

## Comparison note

The model side needs dM/dH, which is a numerical derivative of a computed M(H). At dilution
temperatures the T = 0 minimiser is the right tool (classical statistics is *worse* than assuming
T = 0 in the quantum limit -- see the note in `../../magnetization/ppms_2p5K/PROVENANCE.md`), so
`sv_mvh_curve` applies directly and no finite-temperature path is needed here. That makes this
cheaper to model than the 2.5 K data, not more expensive.

A fine field grid is required, since differentiating a coarsely sampled M(H) will smear exactly the
features being looked for.
