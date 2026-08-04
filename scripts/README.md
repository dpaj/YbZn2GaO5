# scripts/ — what to run, and what it costs

46 entry points at the top level, plus `dev/` (benchmarks and tooling) and `legacy/` (the old
analytical co-fit). This index groups them by **what you are trying to do**, because the
filenames alone do not distinguish a 30-second plot from a 12-hour fit.

Everything writes to `results/`, which is **gitignored** — figures do not survive a clone, so
regenerate rather than hunt for them. Read `CLAUDE.md` first for the scientific state and the
conventions that are easy to get wrong.

## Conventions for the cost column

| tag | meaning |
|---|---|
| **free** | seconds to ~1 min, data or CSV only, no Sunny |
| **cheap** | under ~5 min |
| **medium** | ~5–60 min; run under `julia -t auto` |
| **heavy** | hours; use `nohup`/background, and prefer the DGX |

**Run anything with KPM or realization averaging under `julia -t auto`.** Threading is over
realizations for M(H) and over q for KPM. On Windows, disable Modern Standby before a long
run — a suspended box once logged 37,180 s of wall time for 235 s of compute.

Before any long run, `julia scripts/dev/check_julia_sources.jl scripts src` catches the
`@printf`-with-concatenated-format-string trap, which fails at macro expansion *after* the
compute and has cost 950 s of KPM once already.

---

## Look at the data (no model)

Start here when the question is "what do the measurements actually say?".

| script | cost | what it does |
|---|---|---|
| `plot_mpms3_centering_correction.jl` | free | The MPMS3 centring correction. The 0.42 K file's `Moment (emu)` column is EMPTY; the sample sat 3.06 mm off centre, so the fixed-centre fit under-read by 1.486. Resolves a long-standing scale mystery. |
| `plot_ac_susceptibility_plateau_test.jl` | free | NHMFL SCM1 AC susceptibility to 18 T at 20 mK. Read its header before trusting any AC number: only one of three coils held our sample, and there is **no perpendicular measurement**. |
| `plot_neutron_1d_two_incident_energies.jl` | free | Ei = 3.32 vs 4.65 meV overlay. `YZGO_EI_OVERLAY=full` for the wide view. The cleanest evidence on which cuts carry the 2.08 meV magnet background. |
| `convert_magnetization_dat_to_csv.jl` | free | Quantum Design `.DAT` → the two-column CSV the repo consumes. Handles both MPMS3 (`DC Moment Free Ctr`, `He3 temp`) and DynaCool VSM (`Moment (emu)`). |
| `compare_1d_4p65_3p32_backgrounds.jl` | cheap | Background comparison between the two incident energies. |
| `crystal_field_van_vleck.jl` | cheap | Yb³⁺ level scheme, ground-doublet g-tensor, and χ_VV. Gives χ_VV^zz = 0.0171 ± 0.0007 μB/T — a factor 2.2 below what the M(H) fit wants. |

## The background programme

The fit window is ~70% interpolated background, so this chain is what the neutron uncertainties
rest on. Run in order; Stage 1 for Ei = 4.65 is the production one.

| script | cost | what it does |
|---|---|---|
| `background_stage0_ei332.jl` | cheap | Can Ei = 3.32 determine its own background? **No** — only 8 of 52 energies confirm. Establishes that the answer must come from elsewhere. |
| `background_stage1_ei332.jl` | cheap | Physics-informed background for Ei = 3.32. Validates at 9 T, fails informatively at 14 T. |
| `background_stage1_ei465.jl` | medium | **Production.** Builds the background under 36 defensible choices and writes the envelope table that feeds `sv_load_background_sigma`. Brackets *two families* — smooth (under-subtracts) and min-fed (over-subtracts) — because varying window edges alone would have understated the uncertainty tenfold with every variant wrong identically. |
| `plot_chi2_contributions.jl` | medium | Which energies actually drive the fit. Showed that 35–41% of χ² on the 9 T dispersive cuts comes from the 1.8–2.4 meV band where their modes are *not*. |

## Neutron forward model and comparison

| script | cost | what it does |
|---|---|---|
| `plot_neutron_vs_exp.jl` | medium | The six 1D cuts against background-subtracted experiment. The workhorse. |
| `plot_neutron_parameter_sets.jl` | medium | Same, overplotting several **named** parameter sets — use this to compare candidates in one pass rather than re-running one at a time. |
| `plot_neutron_2d_parameter_sets.jl` | heavy | 2D data-vs-model along the CNCS path. Note the 2D experimental data are **raw, not background-subtracted**, so the ~2 meV feature in them is magnet background. |
| `plot_2d_data_vs_model.jl` | heavy | Repo-native single-parameter-set 2D comparison. |
| `sunny_kpm_1d_disp_grid_2sigmaJ_vs_exp.jl` | heavy | Dispersive component only, deterministic 1D Q-grid. Where q-threading and ground-state reuse were developed. |
| `sunny_plot_kpm_2d.jl` | heavy | Thin 2D driver. |
| `compare_analytical_sunny_dispersion.jl` | medium | Analytical field-polarized dispersion against Sunny. Pins the J1/J2 shell-offset convention. |
| `compare_analytical_vs_sunny_outputs.jl` | free | Post-processes already-saved outputs — no recompute. |
| `export_analytical_2d_model_csv.jl` | cheap | Analytical 2D calculation to long-form CSV. |

## Magnetization

| script | cost | what it does |
|---|---|---|
| `plot_mvh_parameter_sets.jl` | medium | M(H) against experiment for several named parameter sets. |
| `sunny_largecell_mvh_classical.jl` | medium | Classical finite-T M(H,T) for the minimal single-disordered-phase model. **Apply `moment_sign = -1.0`** — Sunny's moment is `-gS`. |
| `plot_largecell_mvh_classical.jl` | free | Plot-only companion to the above. |
| `fit_mvh_only.jl` | medium | Fits M(H) alone and tests whether it needs a Van Vleck term. It does — but 2.2× more than the crystal field allows, and the excess is *not* a digitization artifact (it survives against primary data). |
| `map_mvh_landscape.jl` | heavy | Maps the objective landscape. Established that M(H) constrains **B_sat, not J1 and gzz separately**, and cannot constrain σ_J at all. |
| `check_mvh_convergence.jl` | heavy | Convergence and protocol diagnostics. Source of "12×12×1 is sufficient" and the realization-vs-cell-size cost argument. |

## Optimization

| script | cost | what it does |
|---|---|---|
| `run_neutron_optimization.sh` | heavy | Launches N multi-start Nelder–Mead chains then the follow-up analysis. Defaults N=4, T=32 — **32 threads/chain minimizes latency**, which is what Nelder–Mead needs; more processes maximizes throughput but not convergence speed. |
| `optimize_neutron_neldermead.jl` | heavy | One Nelder–Mead start. Run several concurrently rather than raising iterations. |
| `analyze_neutron_optimum.jl` | medium | Follow-up once the chains finish. |
| `check_background_variance_effect.jl` | heavy | Acceptance test for the background-variance term: per-cut gzz preference with the variance on vs off (~1 h). A **validation of the UQ machinery**, not a parameter refinement. It **failed its own prediction** — the 9 T dispersive cuts do not move — which is how we learned that a χ² *budget* does not tell you what drives a parameter. |
| `plot_background_variance_effect.jl` | free | Plots the above from its CSV — no recompute. Read its header for the refuted-mechanism argument. |
| `scan_gamma_first_parameters.jl` | heavy | Factorized "Γ-first" parameter scan. |
| `plot_gamma_first_scan.jl` | free | Reads the scan CSV only — no recompute. |
| `run_cofit_9T14T.jl` | heavy | The neutron + magnetization co-fit driver. |
| `run_cofit_9T14T_smoke.jl` | cheap | Smoke test for the above. Run this first. |

## Health checks — run these before committing or after pulling

| script | cost | what it does |
|---|---|---|
| `../test/runtests.jl` | cheap | **The one that matters.** Physics-invariant regression tests. Judge by "all green, 0 failures" — do *not* quote a fixed assertion count, it moves with the config count and threading. |
| `dev/check_julia_sources.jl` | free | Pre-flight parse + the `@printf` literal-format trap. Costs nothing; run it before every long job. |
| `check_repo.jl` | free | Repository health. |
| `check_config_loads.jl` | free | Every TOML in `configs/` deep-merges cleanly. |
| `check_neutron_cut_loading.jl` | free | The validation layer loads the same cuts the legacy path does. |
| `check_cross_machine_reproducibility.jl` | medium | The neutron objective against a pinned target (χ²_red = 27.047850). Verified identical across Windows/Linux, Julia 1.12.3/1.12.6, and 32/128 cores to 1.3e−8. **Run this after touching the objective** — it is what makes numbers comparable between the two machines. |
| `check_sunny_form_factor.jl` | cheap | Magnetic form-factor setting. |
| `check_sunny_parameter_mapping.jl` | free | How canonical analytical parameters map into Sunny. |
| `check_analytical_histogram_convergence.jl` | medium | Event-sampling convergence for the analytical 1D histogrammer. |
| `check_analytical_histogram_grid_convergence.jl` | heavy | Deterministic-grid convergence. Relevant to open thread 1 (does the MC Q-sampling mode buy anything over the 81-point grid?). |

## `dev/` — benchmarks and tooling

Machine-specific by nature; see the two-machine table in `CLAUDE.md` before trusting a number
measured on the other box. **KPM is memory-bandwidth bound and the ceiling is machine-specific**:
Windows CPU q-threading saturates near 3×, the DGX reaches 16.6×, and both land near the same
*absolute* throughput. One A100 in Float64 roughly ties well-threaded DGX CPU — **multiple GPUs
are the real lever** (4 concurrent A100s give 3.91×, 97.8% of ideal).

`dev/benchmark_kpm_1d_scaling.jl`, `dev/benchmark_kpm_gpu_vs_cpu.jl`,
`dev/benchmark_sunny_kpm_cpu_gpu_compare.jl`, `dev/benchmark_sunny_kpm_core.jl`,
`dev/benchmark_process_scaling.jl` (+ its `.sh` driver), and `dev/check_julia_sources.jl`.

`sunny_kpm_gpu_fixed_model_summary.jl` at the top level is the GPU fixed-model summary and
belongs to this group in spirit.

## `legacy/` — the old analytical co-fit

Still in use, but **hardcoded `C:\Users\vdp\...` paths make all three fail on Linux.** The fix
is the env-var-with-Windows-fallback pattern already at
`plot_yzgo_2d_data_vs_model_legacy.jl:5910` (`get(ENV, "YZGO_DATA_DIR", raw"C:\...")`). The data
they want live in `YZGO/CNCS_data/`, **outside this repo**, so they must be copied to the DGX
separately.

Note the legacy analytical co-fit uses a **different magnetization normalization** — the Sunny
path divides by `(1+r2)` and the legacy one does not — so `A_M` is not comparable between them.

## Preliminary bridges

`sunny_validate_magnetization_largecell.jl`, `sunny_validate_magnetization_meanfield.jl`,
`sunny_validate_spinwave_kpm.jl` — thin early validation wrappers, superseded by the drivers
above but kept because they are the smallest working examples of each path.

## `extract_neutron_magnetization_features.jl`

Empirical feature extraction for the 1D scans and M(H); feeds `src/feature_extraction.jl`.
