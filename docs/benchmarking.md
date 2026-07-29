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
