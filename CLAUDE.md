# YbZn2GaO5 (YZGO) — orientation

Read this first. It is the handoff document: it records the current scientific
state, the conventions that are easy to get wrong, and the open threads. The
per-topic detail lives in `docs/`.

## The scientific claim

YbZn2GaO5 was published as a **Dirac quantum spin liquid** (Bag et al., PRL 133,
266703 (2024); Wu et al., PRL 135, 046704 (2025)). The argument of this repository
is that the data are instead explained by **quenched disorder** from Zn/Ga site
mixing. Zhao et al. (PRB 113, 014437 (2026)) independently support the disorder
picture, refining site mixing at x = 0.60(5), y = 0.35(5) and explicitly
anticipating a distribution of Lande g factors. The two camps disagree head-on
about whether the crystal-field broadening is site mixing or CEF-phonon coupling.

Status: a minimal single-disordered-phase model reproduces M(H) semi-quantitatively
(rms 0.004 uB against a 1.12 uB signal) and the 1D neutron cuts qualitatively.

## Layout

```text
configs/   TOML run controls, one per script; deep-merged onto a base config
data/      experimental inputs (tracked; data/raw and data/generated are ignored)
docs/      the real documentation - see the reading order below
results/   generated figures and tables — GITIGNORED, all regenerable
scripts/   runnable entry points; scripts/dev holds benchmarks; scripts/legacy the old co-fit
src/       sunny_validation.jl (large), feature_extraction.jl, parameters.jl
test/      runtests.jl — physics-invariant regression tests (see Housekeeping)
```

`../references/` (a sibling of this repo, **not** tracked) holds the five key
papers. A clone will not have them.

Reading order for the current work: `docs/largecell_mvh_classical.md`,
`docs/crystal_field_van_vleck.md`, `docs/benchmarking.md`, then
`docs/companion/` for the older analytical layer.

## Conventions that are easy to get wrong

- **Effective crystal, not the CIF.** `sv_effective_triangle_crystal` builds a
  **P1** one-site triangular net. Because it is P1, Sunny will not generate the
  triangular shells from one representative bond — the J1 and J2 offsets are
  hard-coded in `sv_j1_shell_offsets` / `sv_j2_shell_offsets`. The test suite pins
  them against the analytical J(q) form factors.
- **Isotropic Heisenberg.** J1 and J2 are scalar exchanges. There is no Delta / XXZ
  anisotropy and no single-ion anisotropy anywhere; all anisotropy lives in the
  g-tensor. Note the published fits use an XXZ model with Delta ~ 1.35, so this
  model differs from the literature.
- **`sunny_transverse_gxy = 1.0` is an intensity gauge, not a physical gperp.**
  `ssf_perp` needs nonzero transverse moment components or the inelastic intensity
  vanishes numerically. Never report it as a refined g-factor.
- **Sign.** Sunny's moment is `-g S`, so `sv_m_parallel_uB_per_site` returns the
  calculator sign and `[largecell].moment_sign = -1.0` converts to the
  experimental convention. Always apply it.
- **Disorder.** `sigma_J` is **fractional** (`J*(1 + sigma_J*randn)`), `sigma_gzz`
  is **absolute** and clipped at zero. `sv_apply_disorder!` takes a `realization`
  keyword; `realization = 0` reproduces the original single-realization seed
  bit-for-bit, so published neutron/KPM results are unchanged.
- **Two magnetization normalizations exist.** The Sunny path divides by `(1+r2)`;
  the legacy analytical co-fit does not. `magnetization_global_scale = 0.4189` came
  from the un-normalized fit, so `A_M` is not comparable between them.
- **`[largecell]` in the base config is superseded** except for `moment_sign` and
  the `include_*_disorder` flags. Its `4x4x1` drives only the legacy script.

## Established results — do not re-derive

- **M(H) constrains B_sat, not J1 and gzz separately.** With `A_M` profiled out,
  gzz enters only through `B_sat = S*D_max/(gzz*mu_B)`; its amplitude role is
  absorbed. The ratio is ~2x better determined than either parameter. M(H) wants
  `B_sat ~ 4.0 T` where the by-eye neutron parameters give 5.1 T.
- **M(H) cannot constrain `sigma_J` at all.** Over the whole range 0 to 1 the rms
  moves by 2.7x the reproducibility floor. `sigma_J` must come from the spectra.
  `sigma_gzz` by contrast is tightly constrained. The two observables are therefore
  complementary, which is ideal for a co-optimization.
- **Temperature is irrelevant to M(H) at 0.42 K.** T = 0 and classical 0.42 K give
  rms 0.0182 vs 0.0176. Disorder, not temperature, produces both the rounding of
  the saturation and the extra initial slope, and both are present at T = 0.
- **12x12x1 is sufficient for M(H).** Realization scatter goes as N^-0.38 while
  cost goes as N^1.40, so accuracy is bought with realizations, not cell size
  (3.2x cheaper per unit accuracy). The ground-state texture has a correlation
  length of order one lattice constant, so the box is not limiting.
- **KPM is memory-bandwidth bound.** CPU threading over q saturates at ~3x
  (verified: not BLAS, not construction overhead). Outer-loop CPU parallelism will
  not help either. The GPU port is the only effective lever, and multiple GPUs are
  the way to combine "threading" with GPU.
- **The M(H) linear term is not Van Vleck.** The crystal field gives
  `chi_VV^zz = 0.0171 +- 0.0007 uB/T`; the fit wants 0.0368, a factor 2.2 more.

## Protocol for M(H) work

T = 0 `minimize_energy!` from a field-polarized start with adiabatic field
continuation, 12x12x1, 8-16 **fixed** disorder realizations (common random numbers)
so the objective is deterministic in the parameters. Validate any optimum at
36x36x1 with different seeds. Threading is over realizations, so use `julia -t auto`.

## Gotchas

- Run anything with KPM or realization averaging under **`julia -t auto`**.
- `sv_interp1` returns NaN outside the measured range. The M(H) data span
  0.0225-6.975 T; the objective configs use 0.2-6.8 T for margin, which discards
  usable low-field data (see open threads).
- The KPM kernel FWHM is bounded above by the instrument: the CNCS table *falls*
  with energy transfer (0.155 meV at 0.5 meV to 0.055 at 4 meV). A kernel wider
  than the target cannot be corrected by deposition. Use
  `sv_resolution_kernel_headroom` before a run.
- `[kpm].tol = 0.05` carries ~10% rms error against tol = 0.005. Acceptable for
  looking at a lineshape, not for fitting; 0.01 costs only 1.7x.
- The MC Q-sampling mode at `n_events = 5000` costs ~36 min per 6-cut comparison on
  CPU. The 81-point deterministic grid costs ~35 s. Whether the extra events buy
  anything has never been tested.
- `results/` is gitignored, so figures do not survive a clone. Regenerate them.

## Open threads

1. **Q-convergence for the neutron cuts** — blocks a neutron objective. Does the
   MC mode buy anything over the 81-point grid?
2. **A neutron objective** mirroring `sv_mvh_objective`: the neutron side is a
   mature forward calculator with no residual metric, no realization averaging, no
   threading, and no optimizer. Ground-state reuse and q-threading live in a
   diagnostic script and should move into the library.
3. **Co-optimization** of the spectra and M(H), which is the point of all of the
   above: M(H) pins B_sat and sigma_gzz, the spectra must pin sigma_J.
4. Does the neutron model still need the phenomenological flat component (r2)
   at high disorder? M(H) no longer does.
5. Where does the excess linear M(H) term come from, if not Van Vleck? Suspect the
   normalization of the digitized data, an impurity, or model error.
6. Zero field. The published claim is a zero-field statement and everything here is
   field-polarized. Hard: no polarized reference state, LSWT invalid, and classical
   statistics is a poor approximation at the 0.07 K neutron temperature.

## Two-machine workflow — read this before touching git

This repo is worked on from **two machines that share state only through GitHub**:

- a **Windows desktop** (paths under `C:\Users\vdp\ORNL Dropbox\...`), and
- **`neutrons-dgx01.ornl.gov`**, a *shared* Linux box with A100s, at
  `~/repos/YbZn2GaO5`, driven over VS Code Remote-SSH.

Remote is `github.com/dpaj/YbZn2GaO5`; work goes straight to `main`. There is no
direct link between the two sessions, and **two separate Claude Code sessions run
on the two machines and cannot see each other.** The human is the only channel
between them, and this file is the only shared memory — per-machine Claude memory
directories do not travel.

Rules that follow from that:

- **`git fetch` before reporting or reasoning about repo state.** `git status`
  saying "up to date with 'origin/main'" only compares HEAD to the *locally cached*
  remote ref. A box that has not fetched will report "up to date" while being many
  commits behind. This has already caused confusion once.
- **Commit and push before switching machines or ending a session.** Do not leave
  work stranded on one side.
- **Never assume the local working copy is current** if the other machine may have
  been touched. Ask.
- **Prefer `git pull --ff-only`.** It refuses rather than silently creating a merge
  commit when the two sides have diverged, which is what you want when you cannot
  see the other machine.
- When both sides have work: **pull first, then commit.** Committing first makes the
  box simultaneously ahead and behind for no reason.

### Files that need care across machines

- `CLAUDE.md` and `.claude/settings.json` are **shared project config** — push and
  pull them like code. `.claude/settings.local.json` is machine-specific and is
  ignored by this repo's `.gitignore` (deliberately not by a machine-local global
  ignore, which would not travel).
- **`Manifest.toml`, `envs/sunny-main/Manifest.toml`, `envs/sunny-kpm-gpu/Manifest.toml`
  are tracked and will churn across OSes.** They carry 167-185 `_jll` deps, and the
  CUDA artifacts in the GPU env resolve differently on Linux. `julia_version` is
  pinned at 1.12.3, so a different Julia guarantees a diff. After running
  `Pkg.instantiate()` on a second machine, **do not commit the churn** — use
  `git checkout -- Manifest.toml` — or it ping-pongs forever.
- There is **no `.gitattributes`**, and 111 tracked files carry CRLF in the index.
  Renormalizing is a one-time noisy commit touching all of them, so it must be done
  while *both* machines are clean and pulled immediately on the other side.
- Three `scripts/legacy/*.jl` files have hardcoded `C:\Users\vdp\...` paths and will
  fail on Linux. They are still in use, so the fix is to read an env var with the
  Windows path as fallback, following the pattern already at
  `plot_yzgo_2d_data_vs_model_legacy.jl:5910`
  (`get(ENV, "YZGO_DATA_DIR", raw"C:\...")`). Note the data they want lives in
  `YZGO/CNCS_data/`, **outside this repo**, so it must be copied to the DGX
  separately.
- `../references/` (13 MB of published PDFs) is outside the repo and must stay
  untracked — it is a shared box, and they are copyrighted.

### Which machine for what

KPM is **memory-bandwidth bound**, so CPU threading saturates near 3x and the GPU
is the only real lever. The Windows box has an RTX A2000 (6 GB, FP64 at 1/32 rate);
the DGX has A100s (~5-7x the bandwidth, uncrippled FP64, no OOM ceiling, and
multiple devices for genuine outer-loop parallelism). Heavy KPM belongs on the DGX;
M(H) work is cheap enough anywhere.

## Housekeeping

Run `julia --project=. test/runtests.jl` before committing; it catches convention
regressions. On a machine that has just pulled, also check `julia --version`
against the manifests' 1.12.3 and run `Pkg.instantiate()` first.

**Do not quote a fixed assertion count.** The suite's total is
`39 + (number of files in configs/) + (2 if threaded)`, because the config testset
asserts once per config file and the KPM-threading testset self-skips without
threads. So it reads 53/55 on a 14-config tree and 54/56 on a 15-config one — all
correct. Judge it by "all green, 0 failures", not by a number. Runtime is likewise
machine-dependent (18-21 s on the Windows box, 34 s on the DGX), so it is not a
regression signal either.
