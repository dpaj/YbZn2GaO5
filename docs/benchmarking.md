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
const N_MEASURED_H = 7
const N_MEASURED_K = 7
const N_MEASURED_L = 1

const N_RESOLUTION_H = 5
const N_RESOLUTION_K = 5
const N_RESOLUTION_L = 1
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
