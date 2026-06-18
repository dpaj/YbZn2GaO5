#!/usr/bin/env julia

# Standalone Sunny.jl KPM CPU/GPU benchmark for the YbZn2GaO5
# field-polarized dispersive component.
#
# This script intentionally has no repo/TOML/data dependencies. It can be run
# with ordinary Sunny.jl, where it benchmarks CPU KPM only, or with an
# experimental Sunny branch that provides KernelAbstractions GPU extensions.
#
# Typical use:
#   julia --project=. benchmark_sunny_kpm_cpu_gpu_compare.jl
#
# Optional environment variables:
#   SUNNY_KPM_ENABLE_GPU=1          attempt GPU benchmark
#   SUNNY_KPM_GPU_BACKEND=CUDA      CUDA, Metal, or OpenCL
#   SUNNY_KPM_GPU_BATCHED=1         prefer Sunny.to_device_batched
#   SUNNY_KPM_GPU_PRECISION=Float64 Float64 or Float32
#   SUNNY_KPM_BENCH_LABEL=main      label used in output filenames
#   SUNNY_KPM_OUTPUT_DIR=...        output directory
#   SUNNY_KPM_BASELINE_CSV=...      compare spectra to a previous run

using Sunny
using LinearAlgebra
using Random
using Printf
using Statistics
using Dates

# -----------------------------------------------------------------------------
# Benchmark knobs
# -----------------------------------------------------------------------------

const SEED = 20260618
const LABEL = get(ENV, "SUNNY_KPM_BENCH_LABEL", "sunny")
const OUTPUT_DIR = get(ENV, "SUNNY_KPM_OUTPUT_DIR", "sunny_kpm_cpu_gpu_benchmark_output")
const BASELINE_CSV = get(ENV, "SUNNY_KPM_BASELINE_CSV", "")

const A_ANGSTROM = 3.376
const C_ANGSTROM = 21.96
const SPIN_S = 0.5

# Canonical YZGO best-fit values. Multipliers let this stay useful for stress
# tests without reading project TOMLs.
const GZZ = 3.784584444
const J1_MEV = 0.2299824987
const J2_MEV = 0.008722296371
const SIGMA_GZZ = 0.3859237343
const SIGMA_J = 0.2396873058
const SIGMA_GZZ_MULTIPLIER = 1.0
const SIGMA_J_MULTIPLIER = 1.0

# Nonzero transverse g is needed for the transverse neutron response in ssf_perp.
# This is an intensity gauge for benchmarking, not a fitted physical g_perp.
const SUNNY_TRANSVERSE_GXY = 1.0

const FIELD_T = 9.0
const FIELD_DIRECTION = [0.0, 0.0, 1.0]

# 3x3 magnetic cell repeated to 36x36x1, matching the current YZGO tests.
const DIMS = (3, 3, 1)
const REPEAT_FACTOR = (12, 12, 1)

const INITIALIZE_FIELD_POLARIZED = true
const RELAX_GROUND_STATE = true
const MAXITERS = 5000

const KPM_TOL = 0.05
const KPM_METHOD = :lanczos       # Required by the KA GPU branch.
const KERNEL_FWHM_MEV = 0.08
const ENERGY_MIN_MEV = 0.6
const ENERGY_MAX_MEV = 2.8
const N_ENERGY = 161

# Q grid: 5x5 measured grid crossed with 5x5 momentum-resolution grid.
# L is intentionally a single point for the core benchmark.
const Q_CENTER = [0.0, 1.0, 0.0]
const MEASURED_HALF_WIDTH = [0.05, 0.05, 0.0]
const RESOLUTION_SIGMA = [0.035, 0.035, 0.0]
const RESOLUTION_NSIGMA = 2.0
const N_MEASURED = (7, 7, 1)
const N_RESOLUTION = (5, 5, 1)

const ENABLE_GPU = lowercase(get(ENV, "SUNNY_KPM_ENABLE_GPU", "0")) in ("1", "true", "yes", "on")
const GPU_BACKEND_NAME = get(ENV, "SUNNY_KPM_GPU_BACKEND", "CUDA")
const GPU_BATCHED = lowercase(get(ENV, "SUNNY_KPM_GPU_BATCHED", "1")) in ("1", "true", "yes", "on")
const GPU_PRECISION_NAME = get(ENV, "SUNNY_KPM_GPU_PRECISION", "Float64")

# -----------------------------------------------------------------------------
# Timing / IO helpers
# -----------------------------------------------------------------------------

struct TimerRows
    rows::Vector{NamedTuple}
end
TimerRows() = TimerRows(NamedTuple[])

function timed!(f, timers::TimerRows, implementation::AbstractString, stage::AbstractString;
                q_samples::Int=0, note::AbstractString="")
    t0 = time_ns()
    result = f()
    elapsed = (time_ns() - t0) / 1e9
    push!(timers.rows, (; implementation=String(implementation), stage=String(stage),
                         seconds=elapsed, q_samples=q_samples, note=String(note)))
    @printf("%-10s %-24s %10.3f s", implementation, stage, elapsed)
    q_samples > 0 && @printf("   (%d Q, %.6f s/Q)", q_samples, elapsed / q_samples)
    isempty(note) || @printf("   %s", note)
    println()
    return result
end

function write_profile_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "implementation,stage,seconds,q_samples,seconds_per_q,note")
        for r in rows
            spq = r.q_samples > 0 ? r.seconds / r.q_samples : NaN
            note = replace(r.note, '"' => "'", ',' => ';')
            @printf(io, "%s,%s,%.9g,%d,%.9g,%s\n",
                    r.implementation, r.stage, r.seconds, r.q_samples, spq, note)
        end
    end
    return path
end

function write_spectrum_csv(path::AbstractString, energies, spectrum)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "energy_meV,intensity_avg")
        for i in eachindex(energies)
            @printf(io, "%.10g,%.10g\n", Float64(energies[i]), Float64(spectrum[i]))
        end
    end
    return path
end

function write_qgrid_csv(path::AbstractString, qpts, weights)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "H,K,L,weight")
        for i in eachindex(qpts)
            q = qpts[i]
            @printf(io, "%.10g,%.10g,%.10g,%.10g\n", q[1], q[2], q[3], weights[i])
        end
    end
    return path
end

function read_spectrum_csv(path::AbstractString)
    lines = readlines(path)
    length(lines) >= 2 || error("No data rows in $path")
    energies = Float64[]
    spectrum = Float64[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, ',')
        length(fields) >= 2 || continue
        push!(energies, parse(Float64, fields[1]))
        push!(spectrum, parse(Float64, fields[2]))
    end
    return energies, spectrum
end

function compare_spectra(a, b)
    n = min(length(a), length(b))
    n > 0 || error("Cannot compare empty spectra")
    aa = Float64.(a[1:n])
    bb = Float64.(b[1:n])
    diff = bb .- aa
    denom = sqrt(mean(abs2, aa)) + eps(Float64)
    return (; n=n,
            rmse=sqrt(mean(abs2, diff)),
            relative_rmse=sqrt(mean(abs2, diff)) / denom,
            max_abs_error=maximum(abs.(diff)),
            integrated_a=sum(aa),
            integrated_b=sum(bb),
            relative_integrated_error=(sum(bb)-sum(aa)) / (abs(sum(aa)) + eps(Float64)))
end

function write_comparison_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "comparison,n,rmse,relative_rmse,max_abs_error,integrated_a,integrated_b,relative_integrated_error")
        for r in rows
            @printf(io, "%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                    r.comparison, r.n, r.rmse, r.relative_rmse, r.max_abs_error,
                    r.integrated_a, r.integrated_b, r.relative_integrated_error)
        end
    end
    return path
end

function write_summary_txt(path::AbstractString, profile_rows, spectra, comparisons; extra="")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "Sunny KPM CPU/GPU benchmark")
        println(io, "timestamp = ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
        println(io, "Sunny version = ", try string(pkgversion(Sunny)) catch; "unknown" end)
        println(io, "label = ", LABEL)
        println(io, "enable_gpu = ", ENABLE_GPU)
        println(io, "gpu_backend = ", GPU_BACKEND_NAME)
        println(io, "gpu_batched = ", GPU_BATCHED)
        println(io, "gpu_precision = ", GPU_PRECISION_NAME)
        println(io)
        println(io, "field_T = ", FIELD_T)
        println(io, "dims = ", DIMS)
        println(io, "repeat_factor = ", REPEAT_FACTOR)
        println(io, "system_size = ", (DIMS[1]*REPEAT_FACTOR[1], DIMS[2]*REPEAT_FACTOR[2], DIMS[3]*REPEAT_FACTOR[3]))
        println(io, "sigma_J_used = ", SIGMA_J_MULTIPLIER * SIGMA_J)
        println(io, "sigma_gzz_used = ", SIGMA_GZZ_MULTIPLIER * SIGMA_GZZ)
        println(io, "q_points = ", N_MEASURED[1]*N_MEASURED[2]*N_MEASURED[3]*N_RESOLUTION[1]*N_RESOLUTION[2]*N_RESOLUTION[3])
        println(io, "n_energy = ", N_ENERGY)
        println(io, "kpm_tol = ", KPM_TOL)
        println(io, "kpm_method = ", KPM_METHOD)
        println(io, "kernel_fwhm_meV = ", KERNEL_FWHM_MEV)
        println(io)
        for implementation in unique(r.implementation for r in profile_rows)
            rows = filter(r -> r.implementation == implementation, profile_rows)
            total = sum((r.seconds for r in rows); init=0.0)
            kpm = sum((r.seconds for r in rows if r.stage == "kpm_intensities"); init=0.0)
            println(io, implementation, "_total_seconds = ", total)
            println(io, implementation, "_kpm_seconds = ", kpm)
            println(io, implementation, "_kpm_fraction = ", total > 0 ? kpm / total : NaN)
        end
        println(io)
        for (name, spec) in spectra
            println(io, name, "_spectrum_sum = ", sum(spec))
        end
        println(io)
        for r in comparisons
            println(io, r.comparison, "_relative_rmse = ", r.relative_rmse)
            println(io, r.comparison, "_relative_integrated_error = ", r.relative_integrated_error)
        end
        isempty(extra) || println(io, extra)
    end
    return path
end

# -----------------------------------------------------------------------------
# Model setup
# -----------------------------------------------------------------------------

# Positive representative shell offsets for the triangular P1 cell, matched to
# the analytical YZGO convention:
# Δ1 = 6 - 2[cos(2πH) + cos(2πK) + cos(2π(H+K))]
# Δ2 = 6 - 2[cos(2π(H-K)) + cos(2π(2H+K)) + cos(2π(H+2K))]
j1_shell_offsets() = ([1, 0, 0], [0, 1, 0], [1, 1, 0])
j2_shell_offsets() = ([1, -1, 0], [2, 1, 0], [1, 2, 0])

function field_unit_vector()
    u = Float64.(FIELD_DIRECTION)
    n = norm(u)
    n > 0 || error("FIELD_DIRECTION must be nonzero")
    return u ./ n
end

function set_g_tensor_all_sites!(sys; gzz::Real, gxy::Real)
    for site in eachsite(sys)
        sys.gs[site] = [Float64(gxy) 0.0 0.0;
                        0.0 Float64(gxy) 0.0;
                        0.0 0.0 Float64(gzz)]
    end
    return sys
end

function build_clean_system()
    units = Units(:meV, :angstrom)
    latvecs = lattice_vectors(A_ANGSTROM, A_ANGSTROM, C_ANGSTROM, 90, 90, 120)
    cryst = Crystal(latvecs, [[0.0, 0.0, 0.0]], 1; types=["Yb"])
    sys = System(cryst, [1 => Moment(s=SPIN_S, g=1.0)], :dipole; dims=DIMS, seed=SEED)

    set_g_tensor_all_sites!(sys; gzz=GZZ, gxy=SUNNY_TRANSVERSE_GXY)
    for off in j1_shell_offsets()
        set_exchange!(sys, J1_MEV, Bond(1, 1, off))
    end
    for off in j2_shell_offsets()
        set_exchange!(sys, J2_MEV, Bond(1, 1, off))
    end

    set_field!(sys, collect(field_unit_vector() .* (FIELD_T * units.T)))
    return (; sys, cryst, units)
end

function make_inhomogeneous_with_disorder!(sys)
    sys = to_inhomogeneous(repeat_periodically(sys, REPEAT_FACTOR))
    rng = MersenneTwister(SEED + 7919)

    sigma_J_used = SIGMA_J_MULTIPLIER * SIGMA_J
    for off in j1_shell_offsets()
        for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
            Jij = J1_MEV * (1.0 + sigma_J_used * randn(rng))
            set_exchange_at!(sys, Jij, s1, s2; offset=o)
        end
    end
    for off in j2_shell_offsets()
        for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
            Jij = J2_MEV * (1.0 + sigma_J_used * randn(rng))
            set_exchange_at!(sys, Jij, s1, s2; offset=o)
        end
    end

    sigma_gzz_used = SIGMA_GZZ_MULTIPLIER * SIGMA_GZZ
    for site in eachsite(sys)
        gi = max(0.0, GZZ + sigma_gzz_used * randn(rng))
        sys.gs[site] = [SUNNY_TRANSVERSE_GXY 0.0 0.0;
                        0.0 SUNNY_TRANSVERSE_GXY 0.0;
                        0.0 0.0 gi]
    end
    return sys
end

function initialize_field_polarized!(sys)
    u = field_unit_vector()
    for site in eachsite(sys)
        set_dipole!(sys, u, site)
    end
    return sys
end

function construct_kpm(sys)
    measure = ssf_perp(sys)
    try
        return SpinWaveTheoryKPM(sys; measure=measure, tol=KPM_TOL, method=KPM_METHOD)
    catch err
        # Older/stable Sunny versions may not expose a method keyword. They use
        # the Lanczos KPM path by default in the documented example.
        @warn "SpinWaveTheoryKPM(...; method=...) failed; retrying without method keyword" exception=(err, catch_backtrace())
        return SpinWaveTheoryKPM(sys; measure=measure, tol=KPM_TOL)
    end
end

# -----------------------------------------------------------------------------
# Q grid / intensity helpers
# -----------------------------------------------------------------------------

function linspace_offsets(half_width::Real, n::Int)
    n <= 1 && return [0.0]
    return collect(range(-Float64(half_width), Float64(half_width); length=n))
end

function gaussian_offsets(sigma::Real, nsigma::Real, n::Int)
    if sigma <= 0 || n <= 1 || !isfinite(sigma)
        return [0.0], [1.0]
    end
    xs = collect(range(-Float64(nsigma)*Float64(sigma), Float64(nsigma)*Float64(sigma); length=n))
    ws = exp.(-0.5 .* (xs ./ Float64(sigma)).^2)
    ws ./= sum(ws)
    return xs, ws
end

function build_q_grid()
    mh = linspace_offsets(MEASURED_HALF_WIDTH[1], N_MEASURED[1])
    mk = linspace_offsets(MEASURED_HALF_WIDTH[2], N_MEASURED[2])
    ml = linspace_offsets(MEASURED_HALF_WIDTH[3], N_MEASURED[3])

    rh, wh = gaussian_offsets(RESOLUTION_SIGMA[1], RESOLUTION_NSIGMA, N_RESOLUTION[1])
    rk, wk = gaussian_offsets(RESOLUTION_SIGMA[2], RESOLUTION_NSIGMA, N_RESOLUTION[2])
    rl, wl = gaussian_offsets(RESOLUTION_SIGMA[3], RESOLUTION_NSIGMA, N_RESOLUTION[3])

    qpts = Vector{Vector{Float64}}()
    weights = Float64[]
    measured_weight = 1.0 / (length(mh) * length(mk) * length(ml))

    for dhm in mh, dkm in mk, dlm in ml
        for (i, dhr) in pairs(rh), (j, dkr) in pairs(rk), (k, dlr) in pairs(rl)
            push!(qpts, [Q_CENTER[1] + dhm + dhr,
                         Q_CENTER[2] + dkm + dkr,
                         Q_CENTER[3] + dlm + dlr])
            push!(weights, measured_weight * wh[i] * wk[j] * wl[k])
        end
    end
    weights ./= sum(weights)
    return qpts, weights
end

function intensity_data_matrix(Iraw)
    if Iraw isa AbstractMatrix
        return Matrix(Iraw)
    elseif hasproperty(Iraw, :data)
        return Matrix(getproperty(Iraw, :data))
    elseif hasproperty(Iraw, :intensities)
        return Matrix(getproperty(Iraw, :intensities))
    else
        error("Could not convert intensities output of type $(typeof(Iraw)) to a matrix")
    end
end

function normalize_intensity_matrix(Iraw, nE::Int, nQ::Int)
    I = intensity_data_matrix(Iraw)
    if size(I) == (nE, nQ)
        return I
    elseif size(I) == (nQ, nE)
        return transpose(I) |> Matrix
    else
        error("Unexpected intensities matrix size $(size(I)); expected ($nE,$nQ) or ($nQ,$nE)")
    end
end

# -----------------------------------------------------------------------------
# GPU helpers
# -----------------------------------------------------------------------------

function gpu_precision_type()
    s = lowercase(GPU_PRECISION_NAME)
    s in ("float32", "f32", "single") && return Float32
    s in ("float64", "f64", "double") && return Float64
    error("Unsupported SUNNY_KPM_GPU_PRECISION=$GPU_PRECISION_NAME; use Float32 or Float64")
end

function _require_pkg(uuid::String, name::String)
    # Load optional GPU packages without putting `using CUDA` at top level.
    # This avoids forcing CUDA as a dependency of CPU-only runs, and avoids
    # PowerShell/Julia world-age issues from `@eval using CUDA` inside a function.
    return Base.require(Base.PkgId(Base.UUID(uuid), name))
end

function maybe_make_backend()
    ENABLE_GPU || return nothing

    name = lowercase(GPU_BACKEND_NAME)
    if name == "cuda"
        try
            ka = _require_pkg("63c18a36-062a-441e-b654-da1e3ab1ce7c", "KernelAbstractions")
            cuda = _require_pkg("052768ef-5323-5732-b1bb-66c8b64840ba", "CUDA")

            functional = Base.invokelatest(getproperty(cuda, :functional))
            functional || begin
                @warn "CUDA loaded but CUDA.functional() is false; GPU benchmark disabled"
                return nothing
            end

            # In current CUDA/KernelAbstractions, CUDABackend is exported by CUDA
            # rather than KernelAbstractions. `ka` is still loaded above so the
            # KernelAbstractions extension machinery is available.
            CUDABackend = getproperty(cuda, :CUDABackend)
            return Base.invokelatest(CUDABackend)
        catch err
            @warn "Could not initialize CUDA backend; GPU benchmark disabled" exception=(err, catch_backtrace())
            return nothing
        end
    elseif name == "metal"
        @warn "Metal backend is not wired in this benchmark script yet; use CUDA for the current GPU test"
        return nothing
    elseif name == "opencl"
        @warn "OpenCL backend is not wired in this benchmark script yet; use CUDA for the current GPU test"
        return nothing
    else
        @warn "Unknown SUNNY_KPM_GPU_BACKEND=$GPU_BACKEND_NAME; GPU benchmark disabled"
        return nothing
    end
end

function make_device_kpm(swt, backend)
    precision = gpu_precision_type()
    if GPU_BATCHED && isdefined(Sunny, :to_device_batched)
        return Base.invokelatest(getproperty(Sunny, :to_device_batched), swt, backend; precision=precision)
    elseif isdefined(Sunny, :to_device)
        return Base.invokelatest(getproperty(Sunny, :to_device), swt, backend)
    else
        error("Sunny GPU extension entry points were not found. Expected Sunny.to_device_batched or Sunny.to_device.")
    end
end

# -----------------------------------------------------------------------------
# Main benchmark
# -----------------------------------------------------------------------------

function main()
    timers = TimerRows()
    mkpath(OUTPUT_DIR)

    @printf("Sunny KPM CPU/GPU benchmark\n")
    @printf("---------------------------\n")
    @printf("Sunny version       = %s\n", try string(pkgversion(Sunny)) catch; "unknown" end)
    @printf("label               = %s\n", LABEL)
    @printf("field_T             = %.3f\n", FIELD_T)
    @printf("system_size         = (%d, %d, %d)\n", DIMS[1]*REPEAT_FACTOR[1], DIMS[2]*REPEAT_FACTOR[2], DIMS[3]*REPEAT_FACTOR[3])
    @printf("sigma_J used        = %.9g\n", SIGMA_J_MULTIPLIER * SIGMA_J)
    @printf("sigma_gzz used      = %.9g\n", SIGMA_GZZ_MULTIPLIER * SIGMA_GZZ)
    @printf("grid                = measured %dx%dx%d, resolution %dx%dx%d\n",
            N_MEASURED..., N_RESOLUTION...)
    @printf("enable_gpu          = %s\n", string(ENABLE_GPU))
    @printf("gpu_backend         = %s\n", GPU_BACKEND_NAME)
    @printf("gpu_batched         = %s\n", string(GPU_BATCHED))
    @printf("gpu_precision       = %s\n", GPU_PRECISION_NAME)
    println()

    built = timed!(timers, "shared", "build_clean_system") do
        build_clean_system()
    end

    sys = timed!(timers, "shared", "enlarge_and_disorder") do
        make_inhomogeneous_with_disorder!(built.sys)
    end

    timed!(timers, "shared", "initialize_spins") do
        INITIALIZE_FIELD_POLARIZED ? initialize_field_polarized!(sys) : randomize_spins!(sys)
    end

    if RELAX_GROUND_STATE
        timed!(timers, "shared", "minimize_energy") do
            minimize_energy!(sys; maxiters=MAXITERS)
        end
    end

    swt = timed!(timers, "shared", "construct_kpm") do
        construct_kpm(sys)
    end

    energies = collect(range(ENERGY_MIN_MEV, ENERGY_MAX_MEV; length=N_ENERGY))
    kernel = gaussian(fwhm=KERNEL_FWHM_MEV)

    qpts, weights = timed!(timers, "shared", "build_q_grid") do
        build_q_grid()
    end
    nQ = length(qpts)

    spectra = Dict{String, Vector{Float64}}()

    Icpu = timed!(timers, "cpu", "kpm_intensities"; q_samples=nQ, note="Sunny CPU") do
        intensities(swt, qpts; energies=energies, kernel=kernel, kT=0.0, verbose=false)
    end
    spectra["cpu"] = timed!(timers, "cpu", "average_spectrum"; q_samples=nQ) do
        I = normalize_intensity_matrix(Icpu, length(energies), nQ)
        Vector(I * weights)
    end

    backend = maybe_make_backend()
    if backend !== nothing
        swt_gpu = timed!(timers, "gpu", "transfer_to_device"; note="$(GPU_BACKEND_NAME), $(GPU_PRECISION_NAME), batched=$(GPU_BATCHED)") do
            make_device_kpm(swt, backend)
        end
        Igpu = timed!(timers, "gpu", "kpm_intensities"; q_samples=nQ, note="$(GPU_BACKEND_NAME)") do
            # The KA extension methods are loaded dynamically. Use invokelatest
            # so a freshly-loaded Sunny.intensities(::DeviceKPM, ...) method is
            # callable from this script without hitting Julia world-age issues.
            Base.invokelatest(Sunny.intensities, swt_gpu, qpts; energies=energies, kernel=kernel, kT=0.0, verbose=false)
        end
        spectra["gpu"] = timed!(timers, "gpu", "average_spectrum"; q_samples=nQ) do
            I = normalize_intensity_matrix(Igpu, length(energies), nQ)
            Vector(I * weights)
        end
    end

    comparisons = NamedTuple[]
    if haskey(spectra, "gpu")
        c = compare_spectra(spectra["cpu"], spectra["gpu"])
        push!(comparisons, merge((comparison="cpu_vs_gpu",), c))
    end

    if !isempty(BASELINE_CSV) && isfile(BASELINE_CSV)
        _, baseline = read_spectrum_csv(BASELINE_CSV)
        for (name, spec) in spectra
            c = compare_spectra(baseline, spec)
            push!(comparisons, merge((comparison="baseline_vs_$(name)",), c))
        end
    elseif !isempty(BASELINE_CSV)
        @warn "SUNNY_KPM_BASELINE_CSV was set but file was not found" BASELINE_CSV
    end

    profile_path = joinpath(OUTPUT_DIR, "profile_$(LABEL).csv")
    qgrid_path = joinpath(OUTPUT_DIR, "qgrid_$(LABEL).csv")
    comparison_path = joinpath(OUTPUT_DIR, "comparison_$(LABEL).csv")
    summary_path = joinpath(OUTPUT_DIR, "summary_$(LABEL).txt")

    timed!(timers, "shared", "write_outputs") do
        write_profile_csv(profile_path, timers.rows)
        write_qgrid_csv(qgrid_path, qpts, weights)
        for (name, spec) in spectra
            write_spectrum_csv(joinpath(OUTPUT_DIR, "spectrum_$(LABEL)_$(name).csv"), energies, spec)
        end
        write_comparison_csv(comparison_path, comparisons)
        write_summary_txt(summary_path, timers.rows, spectra, comparisons)
    end

    println()
    println("Saved:")
    println("  ", profile_path)
    for name in sort(collect(keys(spectra)))
        println("  ", joinpath(OUTPUT_DIR, "spectrum_$(LABEL)_$(name).csv"))
    end
    println("  ", qgrid_path)
    println("  ", comparison_path)
    println("  ", summary_path)
    return nothing
end

main()
