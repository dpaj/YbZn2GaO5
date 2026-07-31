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
- **KPM is memory-bandwidth bound, but the ceiling is MACHINE-SPECIFIC.** On the
  Windows box CPU q-threading saturates near 3x; on the DGX it reaches 16.6x at 81
  chunks and is still rising, because a DGX has far more aggregate memory bandwidth.
  Both machines land at a similar *absolute* throughput (~0.05-0.06 s/q), which is
  the bandwidth wall; the speedup *factor* differs only because the serial baselines
  do. An earlier claim here that "the GPU port is the only effective lever" was
  wrong for the DGX: a single A100 in Float64 (0.065 s/q) actually LOSES to
  well-threaded DGX CPU (0.0495 s/q). **Multiple GPUs are the real lever** -- 4
  concurrent A100s give 3.91x aggregate, 97.8% of ideal. Host threading on top of
  one GPU helps ~1.2x at 36x36x1 but is net contention (0.88x) at 12x12x1, so it is
  system-size dependent.
- **Ground-state cost is small; an earlier claim that it was ~76% of runtime was a
  benchmarking error.** That figure was first-call JIT, not compute. Warm it is
  0.010 s at 12x12x1 and 0.18 s at 36x36x1 / 14 T (cold/warm ratio 579x at
  12x12x1), and warm times do scale with system size. There is no Amdahl cap on the
  GPU gain.
- **But `maxiters` is badly mis-set for 9 T.** At 9 T on 36x36x1 the energy is
  converged by ~1000 iterations, yet the minimizer never satisfies its convergence
  test (gradient plateaus near 8e-8 and is not even monotonic), so it burns the cap.
  E/site is identical to 8 decimals for maxiters of 1000, 5000, 20000 and 50000,
  while the time goes 0.89, 4.42, 17.8, 47.1 s. `maxiters = 50000` therefore wastes
  ~46 of 47 seconds at exactly the field where the neutron data lives. 14 T
  converges properly in 134 iterations; 9 T is pathological because the disordered
  system is only partially saturated there. **Judge convergence by E/site, not by
  the returned flag**, and do not raise `maxiters` to chase the flag.
- **`sigma_gzz`, not `sigma_J`, dominates both KPM cost and realization scatter.**
  Confirmed from two orthogonal scans on two machines: at fixed `sigma_J = 0.5`, cost rises
  **2.82x** across `sigma_gzz` 0 -> 1.0 (Windows); at fixed `sigma_gzz = 0.8`, it rises only
  **1.12x** across `sigma_J` 0 -> 1.0 (DGX). The realization floor behaves the same way --
  12.4% shape spread rms at `sigma_J = 0`, rising only ~19% by `sigma_J = 1.0`. The reason is
  a scale comparison: the Zeeman term is 1.7 meV at 9 T and 3.1 meV at 14 T and `sigma_gzz =
  0.8` spreads it by +-0.65 meV at 14 T, whereas the whole exchange bandwidth
  `S*[J(0)-J(q)]` is only ~0.5-1 meV at `J1 = 0.25`, `S = 1/2`. Both of us initially
  predicted 2-3x from `sigma_J`; the magnitude was right and the cause wrong.
  This is physics, not just performance: the g-factor distribution is the dominant disorder
  effect in the spectra, which is the mechanism Zhao et al. anticipated.
- **The neutron realization floor is 12-15%** (shape spread rms, per-cut range 0.080-0.191,
  6 realizations at 36x36x1). It is measured in SPECTRUM units, which is the wrong unit for a
  fit -- the number a fit consumes is the scatter in `chi2_red`, and that is still pending.
  Practical consequence already observed: at 4 realizations the Gamma surface separates
  `chi2_red` 9.4 from 50 trivially but does NOT resolve 9.42 from 9.60.
- **The Gamma cut alone prefers `gzz` = 3.30-3.45 with `sigma_gzz` = 0.80.** Complete 6x6
  surface on `(0,1,0)` at 9 and 14 T, where the exchange cancels identically so J does not
  enter to leading order. The by-eye `gzz = 3.8` scores 21.87 against 9.42, i.e. **2.3x
  worse**, and the published 3.44 sits in the minimum. `sigma_gzz = 0.8` is a genuine
  interior minimum, so the by-eye `sigma_gzz` was right and it is `gzz` that was wrong.
  Caveat: a mean `gzz` of 3.3 with `sigma_gzz = 0.8` is a 24% spread, so the convolved
  lineshape peak is NOT at `gzz*mu_B*B` -- do not compare the fitted mean directly against a
  peak-position estimate.
- **A staged factorization can be globally worse than its own starting point.** `gzz` shifts
  the mode energy at every q, so a `gzz` chosen from Gamma alone cannot see the cost it
  imposes at K and M; a 12x12x1 test landed 73% worse on the six-cut `chi2` while the
  "did Gamma move?" consistency check reported success. Always finish with a joint step, and
  note that an UNWEIGHTED six-cut sum is dominated by whichever cuts fit worst -- the DGX
  measured the same parameter point scoring ~9 on `(0,1,0)` alone and ~216 on all six.
- **The M(H) linear term is not Van Vleck.** The crystal field gives
  `chi_VV^zz = 0.0171 +- 0.0007 uB/T`; the fit wants 0.0368, a factor 2.2 more.

## Protocol for M(H) work

T = 0 `minimize_energy!` from a field-polarized start with adiabatic field
continuation, 12x12x1, 8-16 **fixed** disorder realizations (common random numbers)
so the objective is deterministic in the parameters. Validate any optimum at
36x36x1 with different seeds. Threading is over realizations, so use `julia -t auto`.

## Gotchas

- **Do NOT run KPM under `julia -t auto` on a many-core box.** Intra-process q-threading
  loses efficiency fast, and past a knee it goes *negative*. Measured intra-process
  efficiency (one 81-q spectrum, 36x36x1):

  | chunks/threads | 1 | 2 | 4 | 8 | 16 | 32 | 128 |
  |---|---|---|---|---|---|---|---|
  | Windows, 32 cores | 100% | 92% | 73% | 36% | 22% | 10% | - |
  | DGX, 2x EPYC 7742 | - | - | ~100% | 90% | 73% | - | **5%** |

  On the DGX, **one process with 128 threads is the worst configuration measured** -- 21.96 s
  per spectrum against 11.25 s for the *same code* on 16 threads. The fix is to fan out over
  **processes**, threading only as far as it stays efficient: 32 processes x 4 threads gives
  0.859 spec/s against 0.046 for `-t auto`, an **18.7x** difference. On the Windows box the
  equivalent is ~8 processes x 4 threads, worth ~7x over the 3.4x that 32 threads saturates at.
  The ceiling is process-internal (synchronisation or per-chunk overhead in the q-loop), NOT
  memory bandwidth: per-process throughput stays flat as processes are added, 8x16 and 16x8
  differ by 17% at identical total threads, and NUMA pinning made it *worse*. An earlier claim
  in this file that KPM is "memory-bandwidth bound" was wrong.
  Corollaries: keep total threads at or below the PHYSICAL core count (using SMT siblings
  costs throughput); do not `taskset`-pin (it created stragglers on the nodes holding the
  shared Julia sysimage page cache); use dynamic work assignment, since at 8x16 the batch span
  was 105.6 s against a 70 s mean and one straggler paces everything.
- **`relax_attempts` must stay at 1.** The escalation loop in `sv_kpm_context` breaks on
  `minimize_energy!`'s convergence *flag*, not on a KPM rejection, and at 9 T the flag never
  trips -- so anything above 1 escalates maxiters 1000 -> 4000 -> 16000 unconditionally,
  chasing exactly what this file says to ignore. Measured: **21.8x the relaxation cost for
  5.8e-9 in E/site**, and the spectrum changes by at most 0.12%, ~100x below the realization
  floor. `sv_neutron_curves` defaulted to 3 until 519cf2f and silently cost ~1.55x on every
  neutron evaluation. Raise it only if KPM actually rejects states -- that is a
  *regularization* problem, not a maxiters problem.
- `sv_interp1` returns NaN outside the measured range. The M(H) data span
  0.0225-6.975 T; the objective configs use 0.2-6.8 T for margin, which discards
  usable low-field data (see open threads).
- The KPM kernel FWHM is bounded above by the instrument: the CNCS table *falls*
  with energy transfer (0.155 meV at 0.5 meV to 0.055 at 4 meV). A kernel wider
  than the target cannot be corrected by deposition. Use
  `sv_resolution_kernel_headroom` before a run.
- `[kpm].tol = 0.05` carries ~10% rms error against tol = 0.005. Acceptable for
  looking at a lineshape, not for fitting; 0.01 costs only 1.7x.
- The MC Q-sampling mode at `n_events = 5000` is ~36 min per 6-cut comparison on
  Windows CPU but only ~4 min per cut on one A100, so it is no longer cost-blocked.
  Whether the extra events buy anything over the 81-point grid is still untested.
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

- **Set your git identity on every machine**, or git invents one from
  username@hostname (which is how `vdp@neutrons-dgx01.ornl.gov` reached history):

      git config --global user.name  "Daniel Pajerowski"
      git config --global user.email "daniel@pajerowski.com"

  That address is what most existing commits use and what links to the `dpaj`
  GitHub account. `.mailmap` canonicalizes the historical variants for git tooling,
  but GitHub does not honour `.mailmap` for account linking, so setting this up
  front is the only real fix.
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

Windows box: 32 cores, RTX A2000 (6 GB, FP64 at 1/32 rate). CPU q-threading
saturates near 3x here. DGX: 8x A100-SXM4-40GB, of which **GPUs 4-7 are typically
held by vLLM** (~0.7 GiB free) and 0-3 are usable; CPU q-threading reaches 16.6x.

Measured on the DGX at 36x36x1, 81 q, tol 0.05, kernel 0.05 meV:

| path | s per q |
|---|---|
| CPU serial | 0.774 |
| CPU, 81 chunks | 0.0495 |
| GPU Float64 | 0.0653 |
| GPU Float64, 441-q batch | 0.0488 |
| GPU Float32 | 0.0532 |
| 4 GPUs concurrent | 3.91x aggregate |

So one A100 in Float64 roughly ties well-threaded DGX CPU; **the win is running
several GPUs at once**. Device memory is ~0.85 MiB per q, so the 5000-event MC mode
needs ~4.8 GiB -- impossible on the 6 GB A2000, trivial on a 40 GB A100 (~45,000 q
per card). Heavy KPM and the MC Q-sampling mode belong on the DGX; M(H) work is
cheap enough anywhere. Note the DGX CPU serial figure is not comparable to the
Windows box 0.216 s/q: it was measured with the kpm-gpu fork of Sunny 0.9.0 and on
slower per-core hardware.

## Housekeeping

Run `julia --project=. test/runtests.jl` before committing; it catches convention
regressions. On a machine that has just pulled, also check `julia --version`
against the manifests' 1.12.3 and run `Pkg.instantiate()` first.

**Do not quote a fixed assertion count.** The total is not a constant: the config
testset asserts once per file in `configs/`, and the KPM-threading testset
self-skips without threads (contributing two assertions with `-t auto`, zero
without). So the number moves whenever a config is added or an assertion is written,
and it differs between machines and between threaded and unthreaded runs — observed
values include 53, 54, 55, 56 and 57, all correct. **Judge it by "all green, 0
failures".** Runtime is machine-dependent too (about 20 s on the Windows box, 34 s
on the DGX), so it is not a regression signal either.
