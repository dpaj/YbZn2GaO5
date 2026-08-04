# Sunny KPM benchmarking

This note records how to run the standalone Sunny KPM benchmark scripts in `scripts/dev/`.
These scripts are development tools for timing and CPU/GPU comparison. They do not load the
YZGO experimental data, TOML analysis controls, or plotting workflows.

## Scripts

```text
scripts/dev/benchmark_sunny_kpm_core.jl
scripts/dev/benchmark_sunny_kpm_cpu_gpu_compare.jl
```

`benchmark_sunny_kpm_core.jl` is a compact CPU-only Sunny KPM benchmark.

`benchmark_sunny_kpm_cpu_gpu_compare.jl` is the CPU/GPU comparison benchmark. It runs a CPU
reference calculation first. If GPU mode is enabled and the active Julia environment supports
the Sunny KernelAbstractions/CUDA extension, it then transfers the KPM object to the GPU and
runs the same Q grid on the GPU.

## Julia environments

The benchmark workflow uses separate Julia environments so the production YZGO environment
does not have to switch Sunny versions.

```text
envs/sunny-main/      # registered/main Sunny
envs/sunny-kpm-gpu/   # experimental Sunny branch with KA/CUDA KPM support
```

### Main Sunny environment

From the repository root:

```powershell
julia --project=envs\sunny-main
```

Inside Julia:

```julia
using Pkg
Pkg.status()
```

### GPU branch environment

From the repository root:

```powershell
julia --project=envs\sunny-kpm-gpu
```

Inside Julia:

```julia
using Pkg
Pkg.status()
```

The GPU environment should show `Sunny` installed from the `kpm-gpu` branch of:

```text
https://github.com/rbhirud2005/Sunny.jl
```

It should also include `KernelAbstractions` and `CUDA`.

## CPU baseline from the command line

Use this to get a fresh-process CPU reference. This includes Julia/Sunny/CUDA package loading
outside the timed sections, but timed stages may still include first-call compilation.

```powershell
$env:SUNNY_KPM_BENCH_LABEL="main_cpu"
$env:SUNNY_KPM_ENABLE_GPU="0"

julia --project=envs\sunny-main scripts\dev\benchmark_sunny_kpm_cpu_gpu_compare.jl
```

Typical outputs:

```text
sunny_kpm_cpu_gpu_benchmark_output/profile_main_cpu.csv
sunny_kpm_cpu_gpu_benchmark_output/spectrum_main_cpu_cpu.csv
sunny_kpm_cpu_gpu_benchmark_output/qgrid_main_cpu.csv
sunny_kpm_cpu_gpu_benchmark_output/comparison_main_cpu.csv
sunny_kpm_cpu_gpu_benchmark_output/summary_main_cpu.txt
```

## GPU branch benchmark from the command line

This uses the GPU-branch Sunny environment. The script still runs a CPU reference first, then
runs the GPU calculation if enabled.

```powershell
$env:SUNNY_KPM_BENCH_LABEL="kpm_gpu_branch"
$env:SUNNY_KPM_ENABLE_GPU="1"
$env:SUNNY_KPM_GPU_BACKEND="CUDA"
$env:SUNNY_KPM_GPU_BATCHED="1"
$env:SUNNY_KPM_GPU_PRECISION="Float64"
$env:SUNNY_KPM_BASELINE_CSV="sunny_kpm_cpu_gpu_benchmark_output\spectrum_main_cpu_cpu.csv"

julia --project=envs\sunny-kpm-gpu scripts\dev\benchmark_sunny_kpm_cpu_gpu_compare.jl
```

Typical outputs:

```text
sunny_kpm_cpu_gpu_benchmark_output/profile_kpm_gpu_branch.csv
sunny_kpm_cpu_gpu_benchmark_output/spectrum_kpm_gpu_branch_cpu.csv
sunny_kpm_cpu_gpu_benchmark_output/spectrum_kpm_gpu_branch_gpu.csv
sunny_kpm_cpu_gpu_benchmark_output/qgrid_kpm_gpu_branch.csv
sunny_kpm_cpu_gpu_benchmark_output/comparison_kpm_gpu_branch.csv
sunny_kpm_cpu_gpu_benchmark_output/summary_kpm_gpu_branch.txt
```

## Warm-session timing

For performance comparisons, it is useful to run benchmarks repeatedly in the same Julia
session. A command-line run is reproducible and representative of a fresh process, but the
first timed call can still include method compilation, extension loading, CUDA initialization,
GPU kernel compilation, and device setup. A second `include(...)` in the same Julia session is
usually closer to steady-state throughput.

Start Julia in the GPU environment:

```powershell
julia --project=envs\sunny-kpm-gpu
```

Then run inside Julia:

```julia
ENV["SUNNY_KPM_BENCH_LABEL"] = "kpm_gpu_branch_warmup"
ENV["SUNNY_KPM_ENABLE_GPU"] = "1"
ENV["SUNNY_KPM_GPU_BACKEND"] = "CUDA"
ENV["SUNNY_KPM_GPU_BATCHED"] = "1"
ENV["SUNNY_KPM_GPU_PRECISION"] = "Float64"
ENV["SUNNY_KPM_BASELINE_CSV"] = "sunny_kpm_cpu_gpu_benchmark_output/spectrum_main_cpu_cpu.csv"

include("scripts/dev/benchmark_sunny_kpm_cpu_gpu_compare.jl")
```

Then run a second time in the same session with a new label:

```julia
ENV["SUNNY_KPM_BENCH_LABEL"] = "kpm_gpu_branch_warm"
include("scripts/dev/benchmark_sunny_kpm_cpu_gpu_compare.jl")
```

The second run is the better estimate of warm steady-state throughput for repeated Sunny KPM
calculations in a long-lived Julia process.

## Environment variables

The CPU/GPU comparison script is controlled by environment variables.

| Variable | Typical value | Meaning |
|---|---:|---|
| `SUNNY_KPM_BENCH_LABEL` | `main_cpu`, `kpm_gpu_branch`, `kpm_gpu_branch_warm` | Label used in output filenames |
| `SUNNY_KPM_ENABLE_GPU` | `0` or `1` | Whether to attempt the GPU calculation |
| `SUNNY_KPM_GPU_BACKEND` | `CUDA` | GPU backend to use |
| `SUNNY_KPM_GPU_BATCHED` | `1` | Whether to use the batched GPU device path when available |
| `SUNNY_KPM_GPU_PRECISION` | `Float64` or `Float32` | Device precision |
| `SUNNY_KPM_BASELINE_CSV` | path to `spectrum_main_cpu_cpu.csv` | Optional baseline spectrum for comparison |
| `SUNNY_KPM_TRACK_MEMORY` | `1` or `0` | Include lightweight CPU allocation and CUDA memory snapshots in the profile output; default is `1` |

## Changing the Q-grid size

The default benchmark grid is:

```text
5 × 5 measured grid × 5 × 5 resolution grid × 1 L = 625 Q points
```

To test a larger 1225-Q grid, edit the constants near the top of
`scripts/dev/benchmark_sunny_kpm_cpu_gpu_compare.jl`:

```julia
const N_MEASURED = (7, 7, 1)
const N_RESOLUTION = (5, 5, 1)
```

This gives:

```text
7 × 7 measured grid × 5 × 5 resolution grid × 1 L = 1225 Q points
```

For timing comparisons, change the output label before rerunning:

```julia
ENV["SUNNY_KPM_BENCH_LABEL"] = "kpm_gpu_branch_warm_1225q"
include("scripts/dev/benchmark_sunny_kpm_cpu_gpu_compare.jl")
```

## What to compare

The most important files are:

```text
profile_<label>.csv
comparison_<label>.csv
summary_<label>.txt
```

Useful quantities:

```text
CPU KPM seconds per Q
GPU KPM seconds per Q
CPU/GPU relative RMSE
CPU/GPU relative integrated intensity error
baseline-vs-current relative RMSE
```

## Memory tracking

The profile CSV also includes lightweight memory diagnostics. These are intended for scaling
checks rather than full memory profiling. The key columns are:

```text
cpu_alloc_bytes
cpu_live_after_bytes
gpu_used_after_bytes
gpu_used_delta_bytes
gpu_free_after_bytes
gpu_total_bytes
```

`cpu_alloc_bytes` is the Julia allocation count reported for that timed stage. The CUDA memory
columns are populated only after the CUDA package/backend is loaded. They are most useful for
checking how GPU memory changes when varying Q-grid size, system size, precision, or batched
versus non-batched execution.

To disable memory snapshots for a minimal timing run:

```powershell
$env:SUNNY_KPM_TRACK_MEMORY="0"
```

For the first successful Float64 test, CPU and GPU spectra should agree to numerical precision
before using the GPU branch for scientific diagnostics.

## 1D KPM energy-scan cost model (scripts/dev/benchmark_kpm_1d_scaling.jl)

Run with `julia -t auto --project=. scripts/dev/benchmark_kpm_1d_scaling.jl`.
Measured on the 36x36x1 disordered cell at the by-eye parameters (J1 = 0.25,
sigma_J = 0.5, sigma_gzz = 0.8), 9 T, config tol = 0.05, kernel FWHM = 0.08 meV,
161 energy points over 0-4 meV, 32 cores.

### Cost scalings

| knob | behaviour | lever? |
|---|---|---|
| number of q | strictly linear, 0.207-0.224 s per q | yes, the dominant one |
| `tol` | cost ~ -log10(tol); 0.005 -> 0.05 is only 2.0x | weak |
| kernel FWHM | cost ~ 1/fwhm; 0.08 -> 0.16 is 2.1x | **no, see below** |
| n_energy | 4.5 s at 41 points, 6.1 s at 321 | essentially free |
| ground state | 51-59 s once per field, amortized | no |

`n_energy` being nearly free means the Chebyshev moments dominate and there is no
reason to coarsen the energy axis.

### The kernel FWHM is NOT free speed

Cost really does go as 1/fwhm, but the CNCS resolution table in
`[kpm.energy_resolution].fwhm_table_meV` *decreases* with energy transfer:
0.155 meV at 0.5 meV falling to 0.055 meV at 4 meV. The KPM kernel of 0.08 meV is
therefore **already wider than the instrument above about 2.7 meV**, where the
`subtract_kpm_kernel` quadrature subtraction goes negative and is clamped by
`min_sigma_meV`. Correctness wants a *finer* kernel at the top of the range, which
costs more, not less. Within the [0.5, 3.0] meV fitting window a kernel of about
0.06 meV would be safe, at roughly 1.3x the cost.

### `tol = 0.05` is a bigger approximation than it looks

Against a `tol = 0.005` reference on identical q points:

| tol | speedup | rel rms error | peak error |
|---|---|---|---|
| 0.01 | 1.19x | 4.2% | 31% |
| 0.02 | 1.49x | 5.8% | 39% |
| 0.05 | 2.02x | 10.2% | 55% |
| 0.1 | 2.56x | 14.5% | 74% |

The configured `tol = 0.05` carries about 10% rms and 55% peak error. For eyeballing
a lineshape that is tolerable; for a quantitative fit it is not, and tightening to
0.01 costs only 1.7x.

### CPU threading saturates at ~3x — the calculation is memory-bandwidth bound

Sunny's KPM is serial internally, so q points can be threaded externally (one
`SpinWaveTheoryKPM` per thread; `SpinWaveTheory` clones the system, so this is
safe, and the threaded result is bit-identical to serial).

| chunks | speedup | efficiency |
|---|---|---|
| 2 | 1.75x | 88% |
| 4 | 2.88x | 72% |
| 8 | 2.98x | 37% |
| 16 | 2.95x | 18% |
| 32 | 3.04x | 10% |

Near-ideal at 2 threads, collapsing by 8, flat thereafter. Two candidate
explanations were tested and **excluded**: serial cost is flat from BLAS = 1 to
BLAS = 32 (so KPM is not BLAS-bound and there is no oversubscription), and
`SpinWaveTheoryKPM` construction costs only 0.158 s, contributing ~1.5 s at 32
chunks. What remains is **memory bandwidth**: the Chebyshev recursion streams the
whole 2N x 2N problem once per moment, so a few cores saturate the memory
controller and further cores buy nothing.

Consequences:

- Do not expect more than ~3x from CPU threading, and there is no point giving KPM
  more than ~8 threads.
- Parallelising the outer loop over cuts, fields or realizations instead will *not*
  help, because the ceiling is a machine-level bandwidth limit rather than a
  decomposition problem.
- **This is why the GPU port matters.** GPU memory bandwidth is roughly an order of
  magnitude higher than CPU, which is the likely origin of the 5-8x measured
  elsewhere in this document. For KPM throughput it is the only effective lever.

### CORRECTIONS from DGX measurements (2026-07-29)

Three claims in the section above are wrong or machine-specific. They are left in
place for provenance; read these first.

**Ground-state cost was a benchmarking error.** The 51-59 s figure was first-call
JIT, not compute -- the ground state was timed cold while only the KPM path had been
warmed. Verified independently on the Windows box: 12x12x1 costs 5.79 s cold and
0.010 s warm, a 579x ratio, and warm times do scale with system size (0.010 / 0.102
/ 0.180 s for 144 / 576 / 1296 sites). Warm at 36x36x1 and 14 T is 0.18 s, matching
the DGX 0.20 s. So there is **no Amdahl cap** on the GPU gain, and the Amdahl table
below is void.

**`maxiters` is mis-set for 9 T, and that is the real ground-state cost.** At 9 T on
36x36x1 the energy converges by ~1000 iterations but the minimizer never satisfies
its convergence criterion: the gradient plateaus near 8e-8 and is not monotonic
(7.88e-8 at 20k iterations, 7.98e-8 at 50k). E/site is identical to 8 decimal places
for maxiters 1000, 5000, 20000, 50000 while wall time goes 0.89, 4.42, 17.8, 47.1 s.
So `maxiters = 50000` wastes ~46 of 47 s at one of the two neutron fields. 14 T
converges in 134 iterations. Judge convergence by E/site, not the returned flag.

**CPU threading does not saturate at 3x everywhere.** That is a property of the
Windows box, not of KPM. On the DGX it reaches 16.6x at 81 chunks and is still
rising. Both machines converge to a similar absolute throughput (~0.05-0.06 s/q),
consistent with a shared memory-bandwidth wall; the speedup factors differ because
the serial baselines do (0.216 s/q Windows, 0.774 s/q DGX with the kpm-gpu fork on
slower cores). Consequently the claim that GPU is "the only effective lever" is
false on the DGX, where one A100 in Float64 (0.0653 s/q) loses to well-threaded CPU
(0.0495 s/q). Multiple GPUs are the lever: 4 concurrent A100s give 3.91x, 97.8% of
ideal. Host threading over a single GPU gives ~1.2x at 36x36x1 but costs 0.88x at
12x12x1, so it is system-size dependent.

**Open accuracy question in the kpm-gpu fork.** On A100, GPU Float64 differs from
in-process CPU Float64 by relative rms 2.5e-6, where the RTX A2000 gave 6.7e-10.
Float32 gives 6.1e-6, only 2.4x worse than Float64 where ~100x would be expected.
CPU is bit-reproducible across processes (9.2e-17) and PEDANTIC_MATH leaves the rms
bit-identical, so TF32 is excluded. A precision-limiting step shared by both device
paths is implied. One concrete hypothesis worth testing: the Lanczos spectral bounds
(`eigbounds`, `niters_bounds`) may differ on device, which would shift the Chebyshev
rescaling window and hence the derived moment count M -- a different *approximation*
rather than a rounding difference, which would explain why Float32 and Float64 track
each other instead of differing by their precision ratio. Print `lo`, `hi` and M on
both paths to check. For perspective, 2.5e-6 is four orders of magnitude below the
tol = 0.05 truncation error (~10% rms), so it cannot bias a fit whose residual floor
is ~1e-2 -- but it should be understood before device Float64 is trusted for
anything tighter.

#### Resolved: it is the single-q spectral bounds (DGX, 2026-07-30)

The hypothesis above is confirmed, with a sharper mechanism than "the device computes
slightly different bounds". Both paths call the *same* CPU `Sunny.eigbounds` on the
host LSWT object (`ext/KAExt/KPM/SpinWaveTheoryKPMBatched.jl:95`), so this is not
device arithmetic. What differs is the *scheme*:

- **CPU** (`src/KPM/Lanczos.jl:130-138`) recomputes bounds per chain -- per q and per
  observable -- and fixes M from that chain's own spectral width. Measured at 36x36x1,
  81 q, tol 0.05, kernel 0.05 meV: `de` falls into two clusters (6.2378-6.2629 and
  6.8692-6.8742) giving **three distinct moment counts, 324 / 326 / 358**.
- **Device** computes bounds at ONE representative q (`q_idx = nq / 2`), applies a 4x
  safety factor to size a shared buffer (`max_iters_global` = 1404 for Float64, 1428
  for Float32), then terminates each chain *adaptively*: `niters_eff` min 358, median
  358, max at the cap. So chains the CPU would run at 324/326 get 358, and some run to
  ~4x the CPU count.

Different M is a different *approximation*, not a rounding difference, and the scheme
is shared by both device precisions -- which is exactly why Float32 tracks Float64 at
~2x instead of the ~100x their precision ratio implies.

**The error grows with disorder, as the mechanism predicts.** Single-q bounds get worse
as the spectral range varies more across q. Relative rms against an in-process CPU
Float64 reference at 36x36x1, 81 q:

| sigma_J | Float64 | Float32 |
|---------|---------|---------|
| 0.2397 (canonical) | 2.61e-06 | 6.55e-06 |
| 0.5                | 8.88e-06 | 1.82e-05 |
| 1.0                | 7.25e-05 | 1.25e-04 |

That is a 28x degradation from canonical to sigma_J = 1.0, scaling roughly as
sigma_J^2.3. **Quote 7e-05, not 2.4e-06, as the figure for the fitting regime** -- the
canonical number is optimistic by ~28x where we actually fit. It remains ~140x below
the tol = 0.05 truncation floor, so it is still hygiene rather than a blocker, but the
margin is much smaller than the canonical number suggests and it would not survive a
tighter tol.

Secondary finding: `eigbounds` uses a randomized Lanczos start, so the bounds are not
reproducible call to call -- the device's representative `de` came out 6.7615 in the
Float64 run and 6.8701 in the Float32 run from the same host object. That puts a floor
on cross-path reproducibility independent of everything above.

### How work is launched matters ~19x (DGX, 2026-07-31)

Measured on `neutrons-dgx01`: 2x AMD EPYC 7742, 256 logical CPUs, 128 physical cores,
2 sockets, **8 NUMA nodes of 16 physical cores each**. One 81-q spectrum at 36x36x1,
tol 0.05, kernel 0.05 meV, regularization 1e-5, by-eye set with gzz = 3.35. A file
barrier synchronises the timed phase so JIT and startup are excluded. Unpinned unless
noted.

| config | total threads | s / spectrum | aggregate throughput |
|--------|--------------|--------------|----------------------|
| 1 x 128 | 128 | 21.96 | 0.046 spec/s |
| 1 x 16  | 16  | 11.25 | 0.089 |
| 8 x 16 **pinned** | 128 | 12.72 | 0.455 |
| 8 x 16  | 128 | 11.12 | 0.666 |
| 16 x 8  | 128 | 18.07 | 0.778 |
| **32 x 4** | 128 | 32.39 | **0.859** |
| 32 x 8  | 256 | 51.32 | 0.576 |

Aggregate throughput is the measured quantity and is what the table is for. Note that
**s/spectrum is NOT comparable across rows**: at high process counts it includes
inter-process contention, so it cannot be divided by a serial time to get a
thread-scaling efficiency. An earlier version of this section did exactly that, using
an unmeasured ~130 s single-thread baseline, and reported efficiencies that a later
spot check contradicted -- a single process at 2 threads takes 25.83 s, which cannot
be reconciled with `32 x 4` at 32.39 s if both were pure thread scaling. The
intra-process scaling curve needs its own measurement on an idle box, at fixed
process count and varying threads; until then the mechanism below is inference from
the throughput column alone.

**One process with 128 threads is the worst configuration measured**: 0.046 spec/s
against 0.859 for `32 x 4`, a factor of 18.7. Both of those are single-process-count
comparisons of aggregate throughput, so the factor is sound. It is also slower *per
spectrum* (21.96 s) than a single process at 16 threads (11.25 s) -- a like-for-like
comparison, since both are one uncontended process -- so q-threading does scale
negatively somewhere between 16 and 128 threads. CLAUDE.md's "run anything with KPM or
realization averaging under `julia -t auto`" is correct on a desktop and harmful here.

**The ceiling is process-internal, not the memory system.** Three independent facts
rule out a shared bandwidth limit: per-process throughput stayed flat as processes
were added; `8 x 16` and `16 x 8` differ by 17% at *identical* total thread count;
and pinning made things worse. The limit is synchronisation or per-chunk overhead in
the q-loop. An earlier reading of the single-process 16.6x plateau as "the
memory-bandwidth wall" was wrong.

**Pinning hurts.** `taskset` to one node's 16 physical cores produced stragglers:
`/proc/<pid>/numa_maps` showed only 64.4% of *anonymous* heap on the pinned node
(27.5% on N6), while 76.4% of *file-backed* pages sat on N7 -- the shared Julia
sysimage every process reads. The slow processes were exactly those on nodes 6 and 7
(17.59 s and 14.58 s against 10.45-12.40 s elsewhere). Pinning bought little locality
and removed the ability to migrate away from congested nodes. `numactl` is not
installed, so `--membind` could not be enforced.

**Stragglers pace the batch.** At `8 x 16` pinned the span was 105.6 s against a 70 s
mean, with barrier sync exact to 0.04 s. A fan-out driver should use dynamic work
assignment, not static equal splits.

Rules of thumb for this box:

- Total threads ~= 128 (the physical core count). `32 x 8` = 256 threads was *worse*
  than `32 x 4` = 128, so SMT siblings cost throughput.
- Fan out over processes; thread only to 4-8 per process.
- Do not pin.
- Half the box can stay free for other users at no throughput cost: `32 x 4` uses 128
  of 256 logical CPUs and is the best config measured. Memory there is 60 GB.
- Whether throughput keeps rising as threads-per-process -> 1 (pure task parallelism)
  is UNMEASURED. An earlier extrapolation to ~1.0 spec/s rested on the discredited
  efficiency model above and should not be relied on. `32 x 4` at 0.859 spec/s is the
  best measured point.

### Intra-process thread scaling, measured properly (DGX, 2026-07-31)

Single process, idle box, uncontended, 81 q at 36x36x1, BLAS pinned to 1. Every row is
one uncontended process, so unlike the multi-process table above, per-spectrum times
here ARE comparable and speedup is meaningful.

| threads | chunks | s / spectrum | speedup | efficiency |
|---------|--------|--------------|---------|------------|
| 1   | 1  | 44.85 | 1.00 | 100% |
| 2   | 2  | 23.65 | 1.90 | 95% |
| 4   | 4  | 12.74 | 3.52 | 88% |
| 8   | 8  | 8.56  | 5.24 | 65% |
| 16  | 16 | 7.64  | **5.87** (peak) | 37% |
| 32  | 32 | 8.69  | 5.16 | 16% |
| 64  | 64 | 15.26 | 2.94 | 5% |
| 128 | 81 | 20.18 | 2.22 | 2% |

True single-thread cost is **44.85 s**. Speedup peaks at 5.87x on 16 threads and then
falls. `chunks` is capped at `nq = 81`, so the 128-thread row is really 81 chunks.

**Cause, and a library optimisation.** `_sv_kpm_fill_intensity!` constructs
`SpinWaveTheoryKPM` *inside* the `Threads.@threads` loop, so every `intensities` call
rebuilds one KPM object per chunk. Construction measured at 36x36x1: mean 0.356 s,
range 0.166-0.917 s. At 81 chunks that is 81 constructions against only 44.85/81 =
0.55 s of per-chunk work, so construction dominates -- its cost grows linearly with
chunk count while the work per chunk shrinks, which is the shape above. A harness that
pre-built the per-chunk objects OUTSIDE the timed region reached 16.6x at 81 chunks
where the library peaks at 5.87x, so caching them on the context looks worth up to
~2.8x. (Indicative: that 16.6x was at canonical parameters with Sunny's 1e-8
regularization default.)

Practical: 16 threads minimises latency per evaluation, 2-4 threads maximises
throughput per core, and the process fan-out above beats both.

### Q-sampling and the momentum-resolution quadrature (DGX, 2026-07-31)

gzz = 3.350, (0,1,0) cuts, 4 realizations under common random numbers, regularization
1e-5 pinned, no escalation. `grid` is `analytical_cut_volume_grid`, `MC` is
`analytical_cut_volume_mc` at matched q.

**How the quadrature bug was found.** Grid and MC disagreed by 0.70% at 625 q, and two
samplers of one integral cannot converge to different values. They were convolving
different kernels: `grid_nsigma = 1.5` (truncated then renormalised) against
`resolution_nsigma_clip = 3.0`. Measuring the width `sv_gaussian_grid_axis` actually
realised showed the production setting (n = 3, nsigma = 1.5) was **5.9% too narrow**,
that raising `grid_nsigma` at n = 3 was catastrophic (-55.8% at 3 sigma, since three
nodes at +/-3 sigma put ~98% of the weight at the centre), that adding nodes at
nsigma = 1.5 also degraded it (-19.8% at n = 9), and that 3-node **Gauss-Hermite is
exact** at the same cost. Fixed in `d2166ff` via Golub-Welsch, exact for any n, with the
legacy rule retained as `resolution_quadrature = "truncated_gaussian_grid"` and its
5.9% error pinned by a test.

**Effect of the fix, and where the residual went.**

| sampler | q | pre-fix | post-fix |
|---------|---|---------|----------|
| grid | 81  | 9.007105 | 8.950716 |
| grid | 625 | 9.002427 | 8.893709 |
| MC   | 81  | 8.898857 | 8.898353 |
| MC   | 625 | 8.939591 | 8.939148 |
| **grid625 - MC625** | | **+0.70%** | **-0.49%** |

MC barely moved, since it never used the grid quadrature. The grid moved down and the
sign of the gap **inverted**: pre-fix the grid was the narrower kernel, post-fix it is
exact and MC is marginally narrower. The gap shrank only 28%, so it was tempting to
blame MC's 3 sigma clip -- but widening that clip to 8 sigma moved chi2 by only
**0.024%** (8.939148 -> 8.937000), which **refutes** that explanation. The residual is
MC's own 1/sqrt(N) convergence: at 625 events it still carries ~0.5% error, consistent
with mc81 -> mc625 having moved 0.46%. So the fix is complete, and `grid625` is the
trustworthy value.

**Which axis limits convergence.** Total q is
`(n_measured_h * n_measured_k) * (n_h * n_k)`, so 225 q can be spent either way:

| configuration | q | chi2_red | vs grid625 | compute |
|---------------|---|----------|------------|---------|
| 3x3 measured x 3x3 resolution | 81  | 8.950716 | +0.64% | 108.1 s |
| 3x3 measured x 5x5 resolution | 225 | 8.956555 | +0.71% | 197.4 s |
| **5x5 measured x 3x3 resolution** | **225** | **8.888533** | **-0.058%** | **193.9 s** |
| 5x5 measured x 5x5 resolution | 625 | 8.893709 | -- | 440.5 s |

The **measured** axis is the constraint; extra resolution nodes buy nothing. That is
exactly what Gauss-Hermite predicts, since 3 nodes are already exact for the Gaussian's
second moment.

**Production setting: 3x3 measured x 3x3 Gauss-Hermite resolution = 81 q.**

The convergence ranking is 225 q (0.058% of the 625-q answer) ahead of 81 q (0.64%), and
an earlier version of this section recommended 225 q on that basis. That was convergence
for its own sake. The experimental systematic floor is far above both: roughly 70% of the
`[0.5, 3.0]` meV fit window rests on PCHIP-interpolated background (see CLAUDE.md, "The
data, and how far to trust it"), and the disorder-realization floor is 12-15%. Measured
directly at n = 8: 81 q gives chi2_red 27.0479 against 26.7616 at 225 q, a difference of
0.286 -- well below the ~1.0 chi2_red realization scatter at that n. The two settings are
indistinguishable to the fit, and 81 q costs 1.76x less (466 s against 818 s per six-cut
evaluation). Use 81 q.

`tol` sensitivity closes the same way and more strongly: chi2_red is identical to six
decimals across tol = 0.05 / 0.02 / 0.01, and the DIFFERENCE between two dissimilar
parameter points (gzz 3.40 against 3.80) is identical to six decimals as well
(-0.386301 at both 0.05 and 0.01). Since a fit consumes only differences, tol = 0.05 is
exact for fitting purposes at 1.85x less cost than 0.01.

### Per-evaluation budget for the neutron fit (DGX, 2026-08-04)

One full six-cut `sv_neutron_objective` evaluation, 8 realizations under common random
numbers, 81 q, tol 0.05, regularization 1e-5, maxiters 1000, relax_attempts 1,
36x36x1. Post KPM-operator caching. Uncontended unless noted.

Latency against thread count (225 q, so comparable to the rows below it):

| threads | wall | compute |
|---------|------|---------|
| 8  | 1773 s | 1721 s |
| 16 | 1095 s | 1039 s |
| **32** | **818 s** | 761 s |
| 64 | 875 s | 820 s |

**32 threads minimises latency.** Note this is a DIFFERENT optimum from the throughput
table earlier in this file, which favours many processes at 4 threads. Both are correct:
throughput and latency optimise differently, and which one matters depends on the
optimiser. Nelder-Mead is sequential internally, so **latency sets fit wall time**.

At the production 81 q, 32 threads: **466 s per evaluation** solo. Under 4-way
concurrency (4 chains x 32 threads = 128 physical cores) it degrades only 24%, to 579 s
per chain, with all four chains in flight:

| configuration | wall per evaluation | starts explored |
|---------------|--------------------|-----------------|
| 1 chain  | 466 s | 1 |
| 4 chains | 579 s | 4 |

So **a 4-chain multi-start at 100 evaluations per chain is about 16 hours** -- four
independent starts for 1.24x the cost of one, against 52 h if run sequentially. An
overnight job rather than a weekend one.

`chi2_red` came back identical to six decimals across every thread count (8/16/32/64),
across 4-way concurrency, and across runs days apart. Under CRN the objective is
bit-reproducible, which is the property a simplex method depends on.

### The disorder-realization floor is robust

Shape spread rms across 6 realizations, `[0.5, 3.0]` meV, unit-integral normalised,
sigma_J = 0.5, sigma_gzz = 0.8, mean over the 6 cuts:

| condition | mean floor |
|-----------|------------|
| gzz = 3.35, 81 q  | 0.1362 |
| gzz = 3.80, 81 q  | 0.1373 |
| gzz = 3.80, 625 q | 0.1372 |

So the ~12-15% floor is **independent of gzz and of q count** -- as it should be, since
realization scatter is physical rather than a sampling artefact. Per-cut range is
0.080-0.191, and amplitude spread is far smaller (0.1-4.6%): g-disorder redistributes
weight in energy without changing the integrated moment. Averaging n realizations cuts
this as 1/sqrt(n), but under common random numbers the draw is frozen and the residual
is a shared offset rather than noise the optimiser chases.

### Projected cost of a full 1D comparison

Six cuts (3 qtags x 2 fields), using the measured 3.04x threaded plateau:

| Q-sampling mode | q per cut | total q | serial | threaded |
|---|---|---|---|---|
| deterministic grid | 81 | 486 | 107 s | 35 s |
| MC events (as configured) | 5000 | 30000 | 6599 s | 36 min |

The deterministic grid is affordable as a fit objective; the MC mode as currently
configured is not. **Nobody has yet checked whether 5000 MC events buys anything
over the 81-point grid** — that Q-convergence study is the direct analogue of the
M(H) realization-count study and is the prerequisite for a neutron objective.
