#!/usr/bin/env julia
#
# GPU-vs-CPU KPM benchmark for YbZn2GaO5.
#
# Ground state and KPM are timed SEPARATELY, and KPM construction is always
# outside the timed region, so the reported figure is seconds-per-q for the
# `intensities` call alone. Both are needed: the ground state is CPU-only and
# unported, so lumping them hides the GPU gain behind Amdahl.
#
# Usage from the repo root:
#
#   # CPU-only measurements run in either environment.
#   julia -t auto --project=. scripts/dev/benchmark_kpm_gpu_vs_cpu.jl
#
#   # GPU measurements need the GPU env and a pinned device. Check nvidia-smi
#   # first on a shared box; GPUs held by other jobs will skew every number.
#   CUDA_VISIBLE_DEVICES=2 julia -t 16 --project=envs/sunny-kpm-gpu \
#       scripts/dev/benchmark_kpm_gpu_vs_cpu.jl
#
#   # Alternate controls file, by argument or environment variable.
#   julia --project=. scripts/dev/benchmark_kpm_gpu_vs_cpu.jl configs/other.toml
#
# Results go to stdout as TSV `RESULT<tab>tag<tab>key<tab>value` lines and, if
# [paths].table_subdir is set, to a CSV under results/feature_tables/.
#
# Concurrency measurement (7b): launch one process per device against a shared
# barrier directory so that the KPM phases actually overlap.
#
#   B=$(mktemp -d)
#   for D in 0 1 2 3; do
#     CUDA_VISIBLE_DEVICES=$D SUNNY_KPM_BENCH_MEASUREMENTS=concurrent \
#     SUNNY_KPM_BENCH_BARRIER_DIR=$B SUNNY_KPM_BENCH_NPROC=4 \
#     SUNNY_KPM_BENCH_TAG=dev$D julia -t 8 --project=envs/sunny-kpm-gpu \
#       scripts/dev/benchmark_kpm_gpu_vs_cpu.jl &
#   done; wait

using Printf
using LinearAlgebra
using Statistics
using Random
using Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"));         using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"));   using .SunnyValidation
const SV = SunnyValidation

const DEFAULT_CONTROLS = "configs/sunny_kpm_gpu_benchmark_controls.toml"

# ---------------------------------------------------------------------------
# Controls, with environment variables as optional overrides
# ---------------------------------------------------------------------------

_envstr(key, default) = get(ENV, key, default)

function _env_override!(run::Dict, key::AbstractString, env::AbstractString, parse_as::Symbol)
    raw = get(ENV, env, "")
    isempty(raw) && return run
    run[key] = if parse_as === :float
        parse(Float64, raw)
    elseif parse_as === :int
        parse(Int, raw)
    elseif parse_as === :bool
        lowercase(raw) in ("1", "true", "yes", "on")
    elseif parse_as === :intlist
        [parse(Int, strip(s)) for s in split(raw, ",") if !isempty(strip(s))]
    elseif parse_as === :strlist
        [String(strip(s)) for s in split(raw, ",") if !isempty(strip(s))]
    else
        raw
    end
    return run
end

_get(d::Dict, k, default) = haskey(d, k) ? d[k] : default
_tuple3(v) = (Int(v[1]), Int(v[2]), Int(v[3]))

function load_benchmark_controls()
    loaded = SV.sv_load_diagnostic_controls(REPO_ROOT, DEFAULT_CONTROLS;
                                            env_var="SUNNY_KPM_BENCHMARK_CONTROLS")
    (; diag, controls, diag_path, base_path) = loaded
    run = get!(diag, "run", Dict{String,Any}())
    gpu = get!(diag, "gpu", Dict{String,Any}())
    conc = get!(diag, "concurrent", Dict{String,Any}())

    _env_override!(run, "measurements", "SUNNY_KPM_BENCH_MEASUREMENTS", :strlist)
    _env_override!(run, "field_T", "SUNNY_KPM_BENCH_FIELD_T", :float)
    _env_override!(run, "nq_side", "SUNNY_KPM_BENCH_NQ_SIDE", :int)
    _env_override!(run, "cpu_chunks", "SUNNY_KPM_BENCH_CPU_CHUNKS", :intlist)
    _env_override!(run, "host_chunks", "SUNNY_KPM_BENCH_HOST_CHUNKS", :intlist)
    _env_override!(run, "realization", "SUNNY_KPM_BENCH_REALIZATION", :int)
    # Canonical-parameter overrides, e.g. "sigma_J=0.5,J1_meV=0.25". Useful for
    # re-measuring accuracy at the strong disorder we actually fit at, where the
    # device's single-q spectral bounds are expected to be a worse approximation.
    let raw = get(ENV, "SUNNY_KPM_BENCH_PARAM_OVERRIDES", "")
        if !isempty(raw)
            d = get!(run, "param_overrides", Dict{String,Any}())
            for kv in split(raw, ",")
                isempty(strip(kv)) && continue
                k, v = split(strip(kv), "=")
                d[String(strip(k))] = parse(Float64, strip(v))
            end
        end
    end
    _env_override!(gpu, "precisions", "SUNNY_KPM_BENCH_PRECISIONS", :strlist)
    _env_override!(gpu, "pedantic_math", "SUNNY_KPM_BENCH_PEDANTIC", :bool)
    _env_override!(conc, "barrier_dir", "SUNNY_KPM_BENCH_BARRIER_DIR", :str)
    _env_override!(conc, "nproc", "SUNNY_KPM_BENCH_NPROC", :int)
    _env_override!(conc, "tag", "SUNNY_KPM_BENCH_TAG", :str)

    return (; diag, controls, run, gpu, conc, diag_path, base_path)
end

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

const ROWS = Vector{NamedTuple}()

function emit(tag, key, value)
    println("RESULT\t", tag, "\t", key, "\t", value)
    push!(ROWS, (; tag = String(tag), key = String(key), value = string(value)))
    flush(stdout)
end

# ---------------------------------------------------------------------------
# System construction, mirroring scripts/sunny_kpm_gpu_fixed_model_summary.jl
# ---------------------------------------------------------------------------

function build_ground_state(params, controls::Dict, run::Dict, field_T::Real;
                            maxiters::Union{Nothing,Integer}=nothing)
    sizectl = SV.sv_system_size_controls(controls, "kpm")
    kc = controls["kpm"]
    iters = maxiters === nothing ? Int(kc["maxiters"]) : Int(maxiters)
    realization = Int(_get(run, "realization", 0))

    t_build = @elapsed begin
        base = SV.sv_build_effective_sunny_system(params, controls;
                                                 component = :dispersive,
                                                 dims = sizectl.dims, field_T = field_T)
        sys = base.sys
        sys = sizectl.repeat_factor != (1, 1, 1) ?
              to_inhomogeneous(repeat_periodically(sys, sizectl.repeat_factor)) :
              to_inhomogeneous(sys)
        SV.sv_apply_disorder!(sys, params, controls; component = :dispersive,
                              include_exchange = _get(kc, "include_exchange_disorder", true),
                              include_gzz = _get(kc, "include_gzz_disorder", true),
                              realization = realization)
    end

    u = SV.sv_field_direction(controls)
    t_init = @elapsed for site in eachsite(sys)
        set_dipole!(sys, u, site)
    end
    t_min = @elapsed minimize_energy!(sys; maxiters = iters)

    return sys, (; t_build, t_init, t_min, total = t_build + t_init + t_min,
                   maxiters = iters, e_per_site = energy_per_site(sys),
                   nsites = length(collect(eachsite(sys))),
                   system_size = sizectl.system_size)
end

# `regularization` must be forwarded: Sunny defaults to 1e-8, but disorder puts
# magnon modes near zero energy, the BdG matrix loses positive-definiteness and KPM
# aborts with "Not an energy-minimum". 36x36x1 happens to survive 1e-8; smaller
# cells and larger sigma_J do not. See [kpm].regularization in the base config.
function make_kpm(sys, controls::Dict)
    kc = controls["kpm"]
    reg = _get(kc, "regularization", nothing)
    kwargs = (; measure = SV.sv_sunny_measure(sys, controls),
                tol = Float64(kc["tol"]),
                method = Symbol(_get(kc, "method", "lanczos")))
    return reg === nothing ? SpinWaveTheoryKPM(sys; kwargs...) :
           SpinWaveTheoryKPM(sys; kwargs..., regularization = Float64(reg))
end

function build_qs(run::Dict)
    n = Int(_get(run, "nq_side", 9))
    ctr = [Float64(x) for x in _get(run, "q_center", [0.5, 0.0, 0.0])]
    half = Float64(_get(run, "grid_nsigma", 2.0)) * Float64(_get(run, "sigma_rlu", 0.035))
    n == 1 && return [copy(ctr)]
    hs = range(ctr[1] - half, ctr[1] + half; length = n)
    ks = range(ctr[2] - half, ctr[2] + half; length = n)
    return [[h, k, ctr[3]] for h in hs for k in ks]
end

function energy_grid(controls::Dict)
    kc = controls["kpm"]
    return collect(range(Float64(kc["energy_min_meV"]), Float64(kc["energy_max_meV"]);
                         length = Int(kc["n_energy"]))),
           gaussian(fwhm = Float64(kc["kernel_fwhm_meV"]))
end

function call_intensities(obj, qs, energies, kernel)
    try
        return Base.invokelatest(Sunny.intensities, obj, qs;
                                 energies = energies, kernel = kernel,
                                 kT = 0.0, verbose = false)
    catch
        return Base.invokelatest(Sunny.intensities, obj, qs;
                                 energies = energies, kernel = kernel)
    end
end

function as_array(res)
    d = hasproperty(res, :data) ? res.data : SV.sv_try_extract_sunny_intensity(res)
    return Array(d)          # (nE, nq); q is the second axis
end

function relrms(a, b)
    x = vec(Float64.(a)); y = vec(Float64.(b))
    length(x) == length(y) || error("relrms: size mismatch")
    den = sqrt(mean(abs2, y))
    return den == 0 ? NaN : sqrt(mean(abs2, x .- y)) / den
end

# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------

function measure_ground_state(tag, params, controls, run)
    fields = [Float64(x) for x in _get(run, "ground_state_fields_T", [Float64(_get(run, "field_T", 14.0))])]
    iters = [Int(x) for x in _get(run, "ground_state_maxiters", [Int(controls["kpm"]["maxiters"])])]

    # One throwaway build so the reported numbers are compute, not first-call JIT.
    # Without this the cold figure is dominated by compilation and does not scale
    # with system size at all.
    _, cold = build_ground_state(params, controls, run, first(fields); maxiters = 50)
    emit(tag, "gs_cold_total_s", @sprintf("%.3f", cold.total))
    emit(tag, "nsites", cold.nsites)
    emit(tag, "system_size", string(cold.system_size))

    for B in fields, it in iters
        _, gs = build_ground_state(params, controls, run, B; maxiters = it)
        emit(tag, @sprintf("gs_warm_field%.1f_maxiters%d_build_s", B, it), @sprintf("%.3f", gs.t_build))
        emit(tag, @sprintf("gs_warm_field%.1f_maxiters%d_minimize_s", B, it), @sprintf("%.3f", gs.t_min))
        emit(tag, @sprintf("gs_warm_field%.1f_maxiters%d_total_s", B, it), @sprintf("%.3f", gs.total))
        # Judge convergence by E/site, not by the minimizer's returned flag.
        emit(tag, @sprintf("gs_warm_field%.1f_maxiters%d_E_per_site", B, it), @sprintf("%.10f", gs.e_per_site))
    end
end

function measure_cpu_scan(tag, sys, controls, run, qs, energies, kernel)
    nq = length(qs)
    blas0 = BLAS.get_num_threads()
    emit(tag, "blas_threads_default", blas0)

    warm = make_kpm(sys, controls)
    call_intensities(warm, qs[1:1], energies, kernel)

    tt = @elapsed as_array(call_intensities(warm, qs, energies, kernel))
    emit(tag, "cpu_serial_blasdefault_s_per_q", @sprintf("%.4f", tt / nq))

    BLAS.set_num_threads(Int(_get(run, "cpu_scan_blas_threads", 1)))
    emit(tag, "blas_threads_scan", BLAS.get_num_threads())

    ref = nothing
    t1 = NaN
    for nc in [Int(x) for x in _get(run, "cpu_chunks", [1])]
        if nc > Threads.nthreads()
            emit(tag, "skip_cpu_chunks_gt_threads", nc)
            continue
        end
        parts = [qs[i:nc:end] for i in 1:nc]
        outs = Vector{Any}(undef, nc)
        swts = [make_kpm(sys, controls) for _ in 1:nc]
        for s in swts
            call_intensities(s, qs[1:1], energies, kernel)
        end
        t = @elapsed begin
            tasks = map(1:nc) do i
                Threads.@spawn begin
                    outs[i] = as_array(call_intensities(swts[i], parts[i], energies, kernel))
                end
            end
            foreach(wait, tasks)
        end
        emit(tag, "cpu_chunks$(nc)_s_per_q", @sprintf("%.4f", t / nq))
        if nc == 1
            t1 = t; ref = outs[1]
        else
            emit(tag, "cpu_chunks$(nc)_speedup_vs_serial", @sprintf("%.3f", t1 / t))
            if ref !== nothing
                recon = similar(ref)
                for i in 1:nc
                    recon[:, i:nc:end] = outs[i]
                end
                # q points are independent and SpinWaveTheory clones the system,
                # so threading must be bit-identical, not merely close.
                emit(tag, "cpu_chunks$(nc)_bitidentical", recon == ref)
            end
        end
        swts = nothing; GC.gc()
    end
    BLAS.set_num_threads(blas0)
    return ref
end

function measure_gpu(tag, sys, controls, run, gpu, conc, qs, energies, kernel;
                     mode::Symbol = :single, do_host_chunks::Bool = false)
    nq = length(qs)
    CUDAmod = Base.require(Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA"))
    Base.require(Base.PkgId(Base.UUID("63c18a36-062a-441e-b654-da1e3ab1ce7c"), "KernelAbstractions"))
    Base.invokelatest(getproperty(CUDAmod, :functional)) ||
        error("CUDA is not functional; run the CPU measurements under --project=. instead")

    if _get(gpu, "pedantic_math", false)
        Base.invokelatest(getproperty(CUDAmod, :math_mode!), getproperty(CUDAmod, :PEDANTIC_MATH))
    end
    emit(tag, "cuda_device", string(Base.invokelatest(getproperty(CUDAmod, :name),
                                                      Base.invokelatest(getproperty(CUDAmod, :device)))))
    emit(tag, "cuda_visible_devices", _envstr("CUDA_VISIBLE_DEVICES", "<unset>"))
    # CUDA is loaded through Base.require, so the MathMode enum's show method is
    # not in this world age; stringify via invokelatest.
    emit(tag, "cuda_math_mode",
         try
             Base.invokelatest(string, Base.invokelatest(getproperty(CUDAmod, :math_mode)))
         catch
             "unavailable"
         end)

    backend = Base.invokelatest(getproperty(CUDAmod, :CUDABackend))
    isdefined(Sunny, :to_device_batched) ||
        error("Sunny.to_device_batched not found; this needs the kpm-gpu branch (envs/sunny-kpm-gpu)")

    total_mem = Base.invokelatest(getproperty(CUDAmod, :total_memory))
    avail() = Base.invokelatest(getproperty(CUDAmod, :available_memory))
    used_mem() = Base.invokelatest(getproperty(CUDAmod, :used_memory))
    cached_mem() = Base.invokelatest(getproperty(CUDAmod, :cached_memory))
    sync() = Base.invokelatest(getproperty(CUDAmod, :synchronize))
    reclaim() = Base.invokelatest(getproperty(CUDAmod, :reclaim))

    # Sampling the pool needs a spare thread; with -t 1 the peak reads as zero.
    function timed_gpu(prec::DataType, nchunks::Int, label::String)
        parts = [qs[i:nchunks:end] for i in 1:nchunks]
        outs = Vector{Any}(undef, nchunks)
        devs = [Base.invokelatest(Sunny.to_device_batched, make_kpm(sys, controls),
                                  backend; precision = prec) for _ in 1:nchunks]
        for d in devs
            call_intensities(d, qs[1:1], energies, kernel)
        end
        sync(); GC.gc(); reclaim()

        stop = Ref(false); pk_drv = Ref(0); pk_live = Ref(0); pk_res = Ref(0)
        sampler = Threads.@spawn begin
            while !stop[]
                pk_drv[] = max(pk_drv[], total_mem - avail())
                u = used_mem(); c = cached_mem()
                pk_live[] = max(pk_live[], u)
                pk_res[] = max(pk_res[], u + c)
                sleep(0.02)
            end
        end
        t = @elapsed begin
            if nchunks == 1
                outs[1] = as_array(call_intensities(devs[1], parts[1], energies, kernel))
            else
                tasks = map(1:nchunks) do i
                    Threads.@spawn begin
                        outs[i] = as_array(call_intensities(devs[i], parts[i], energies, kernel))
                    end
                end
                foreach(wait, tasks)
            end
            sync()
        end
        stop[] = true; wait(sampler)
        emit(tag, "$(label)_s_per_q", @sprintf("%.4f", t / nq))
        emit(tag, "$(label)_peak_pool_live_MiB", @sprintf("%.1f", pk_live[] / 2^20))
        emit(tag, "$(label)_peak_driver_MiB", @sprintf("%.1f", pk_drv[] / 2^20))
        return t, outs
    end

    if mode === :concurrent
        # Rendezvous so that every process's KPM phase overlaps; without it the
        # ground-state builds stagger and the devices never contend.
        bdir = String(_get(conc, "barrier_dir", ""))
        isempty(bdir) && error("[concurrent].barrier_dir must be set for the concurrent measurement")
        nproc = Int(_get(conc, "nproc", 1))
        mkpath(bdir)
        touch(joinpath(bdir, "ready_$(tag)"))
        t_wait = @elapsed while count(startswith("ready_"), readdir(bdir)) < nproc
            sleep(0.05)
        end
        emit(tag, "barrier_wait_s", @sprintf("%.3f", t_wait))
        timed_gpu(Float64, 1, "gpu_f64_concurrent")
        return nothing
    end

    precs = [p == "Float32" ? Float32 : Float64 for p in _get(gpu, "precisions", ["Float64"])]

    t_cold = @elapsed Base.invokelatest(Sunny.to_device_batched, make_kpm(sys, controls),
                                        backend; precision = Float64)
    emit(tag, "gpu_transfer_cold_s", @sprintf("%.3f", t_cold))
    t_warm = @elapsed Base.invokelatest(Sunny.to_device_batched, make_kpm(sys, controls),
                                        backend; precision = Float64)
    emit(tag, "gpu_transfer_warm_s", @sprintf("%.3f", t_warm))
    GC.gc(); reclaim()

    cpu_ref = nothing; t_cpu = NaN
    if _get(gpu, "inprocess_cpu_reference", true)
        swt = make_kpm(sys, controls)
        call_intensities(swt, qs[1:1], energies, kernel)
        t_cpu = @elapsed cpu_ref = as_array(call_intensities(swt, qs, energies, kernel))
        emit(tag, "cpu_inprocess_serial_s_per_q", @sprintf("%.4f", t_cpu / nq))
        swt = nothing; GC.gc()
    end

    times = Dict{DataType,Float64}()
    for prec in precs
        lbl = prec === Float32 ? "gpu_f32" : "gpu_f64"
        t, outs = timed_gpu(prec, 1, lbl)
        times[prec] = t
        !isnan(t_cpu) && emit(tag, "$(lbl)_speedup_vs_cpu_serial_inprocess", @sprintf("%.2f", t_cpu / t))
        if cpu_ref !== nothing
            emit(tag, "$(lbl)_relrms_vs_cpu_inprocess", @sprintf("%.3e", relrms(outs[1], cpu_ref)))
        end
    end
    haskey(times, Float32) && haskey(times, Float64) &&
        emit(tag, "gpu_f32_speedup_vs_f64", @sprintf("%.3f", times[Float64] / times[Float32]))

    # Host-side threading over q, layered on the device batching. Done here so
    # timed_gpu is in scope; each host thread needs its own device copy, so the
    # peak-memory figures for these rows are expected to grow with nchunks.
    if do_host_chunks && haskey(times, Float64)
        for nc in [Int(x) for x in _get(run, "host_chunks", Int[])]
            if nc > Threads.nthreads()
                emit(tag, "skip_host_chunks_gt_threads", nc)
                continue
            end
            t, _ = timed_gpu(Float64, nc, "gpu_f64_hostchunks$(nc)")
            emit(tag, "gpu_f64_hostchunks$(nc)_speedup_vs_1",
                 @sprintf("%.3f", times[Float64] / t))
        end
    end

    return times
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    loaded = load_benchmark_controls()
    (; controls, run, gpu, conc, diag_path, base_path) = loaded

    tag = let t = String(_get(conc, "tag", ""))
        isempty(t) ? "bench" : t
    end
    measurements = [String(m) for m in _get(run, "measurements", ["cpu_scan"])]

    emit(tag, "controls", relpath(diag_path, REPO_ROOT))
    emit(tag, "base_controls", relpath(base_path, REPO_ROOT))
    emit(tag, "julia_version", string(VERSION))
    emit(tag, "julia_threads", Threads.nthreads())
    emit(tag, "measurements", join(measurements, ","))

    params = SV.sv_load_params(REPO_ROOT, controls).params
    params, applied_overrides = SV.sv_apply_param_overrides(params, run)
    for line in applied_overrides
        emit(tag, "param_override", line)
    end
    qs = build_qs(run)
    energies, kernel = energy_grid(controls)
    kc = controls["kpm"]
    emit(tag, "nq", length(qs))
    emit(tag, "n_energy", Int(kc["n_energy"]))
    emit(tag, "kernel_fwhm_meV", Float64(kc["kernel_fwhm_meV"]))
    emit(tag, "tol", Float64(kc["tol"]))
    emit(tag, "maxiters", Int(kc["maxiters"]))
    emit(tag, "realization", Int(_get(run, "realization", 0)))

    "ground_state" in measurements && measure_ground_state(tag, params, controls, run)

    needs_sys = any(m -> m in ("cpu_scan", "gpu", "host_chunks", "concurrent"), measurements)
    if needs_sys
        field_T = Float64(_get(run, "field_T", 14.0))
        emit(tag, "field_T", field_T)
        sys, gs = build_ground_state(params, controls, run, field_T)
        emit(tag, "gs_for_kpm_total_s", @sprintf("%.3f", gs.total))
        emit(tag, "gs_for_kpm_E_per_site", @sprintf("%.10f", gs.e_per_site))

        "cpu_scan" in measurements &&
            measure_cpu_scan(tag, sys, controls, run, qs, energies, kernel)

        if "concurrent" in measurements
            measure_gpu(tag, sys, controls, run, gpu, conc, qs, energies, kernel; mode = :concurrent)
        elseif "gpu" in measurements || "host_chunks" in measurements
            measure_gpu(tag, sys, controls, run, gpu, conc, qs, energies, kernel;
                        do_host_chunks = ("host_chunks" in measurements))
        end
    end

    subdir = String(_get(get(controls, "paths", Dict{String,Any}()), "table_subdir", ""))
    if !isempty(subdir)
        outdir = SV.sv_repo_path(REPO_ROOT, joinpath("results", "feature_tables", subdir))
        mkpath(outdir)
        path = joinpath(outdir, "kpm_gpu_vs_cpu_$(tag).csv")
        SV.sv_write_rows_csv(path, ROWS)
        println("Wrote ", relpath(path, REPO_ROOT))
    end
    println("### BENCHMARK_DONE ", tag)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
