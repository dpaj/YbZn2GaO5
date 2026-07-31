#!/usr/bin/env julia
#
# Multi-process throughput scaling for the KPM neutron path.
#
# Hypothesis under test: q-threading inside ONE process saturates near 16x because
# that is one NUMA node's worth of cores (this box is 2x EPYC 7742 = 8 NUMA nodes
# x 16 physical cores). If so, the lever is N processes each pinned to one node,
# and aggregate throughput should approach n_nodes x per-node throughput -- NOT
# because more bandwidth appears, but because each process gets LOCAL bandwidth
# instead of paying remote-access penalties across sockets.
#
# A file barrier makes every process enter its timed loop together, so the wall
# span is comparable across processes and startup/JIT is excluded.
#
#   MP_NSPEC    spectra per process in the timed loop (default 6)
#   MP_NQ_SIDE  grid side: 3 -> 81 q (production), 5 -> 625 q
#   MP_BARRIER  rendezvous directory
#   MP_NPROC    number of processes expected at the barrier
#   MP_TAG      label

using Printf, LinearAlgebra, Statistics, Random

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"));        using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"));  using .SunnyValidation
const SV = SunnyValidation

BLAS.set_num_threads(1)      # library threads over q; BLAS threads only contend

getenv(k, d) = get(ENV, k, d)
const NSPEC  = parse(Int, getenv("MP_NSPEC", "6"))
const NQSIDE = parse(Int, getenv("MP_NQ_SIDE", "3"))
const TAG    = getenv("MP_TAG", "p0")
const NPROC  = parse(Int, getenv("MP_NPROC", "1"))

say(args...) = (println(args...); flush(stdout))

controls = SV.sv_load_controls(REPO_ROOT)
kc = controls["kpm"]
kc["dims"] = [3, 3, 1]; kc["system_size"] = [36, 36, 1]
kc["method"] = "lanczos"; kc["tol"] = 0.05; kc["maxiters"] = 1000
kc["regularization"] = 1e-6
kc["energy_min_meV"] = 0.1; kc["energy_max_meV"] = 4.2
kc["n_energy"] = 241; kc["kernel_fwhm_meV"] = 0.05
kc["thread_max_chunks"] = Threads.nthreads()
kc["thread_min_q_per_chunk"] = 1
kc["thread_min_q"] = 1

eh = get!(kc, "experimental_histogram", Dict{String,Any}())
eh["enabled"] = true; eh["mode"] = "analytical_cut_volume_grid"
eh["measured_grid_mode"] = "uniform_grid"
eh["n_measured_h"] = NQSIDE; eh["n_measured_k"] = NQSIDE; eh["n_measured_l"] = 1
qa = get!(kc, "q_averaging", Dict{String,Any}())
qa["enabled"] = true; qa["mode"] = "gaussian_grid"
qa["n_h"] = NQSIDE; qa["n_k"] = NQSIDE; qa["n_l"] = 1

p0 = SV.sv_load_params(REPO_ROOT, controls).params
# Representative production point. Cost is sigma_gzz-dominated, so sigma_gzz = 0.8
# is what makes this a fair proxy; gzz is set near the Gamma-first minimum.
ov = Dict{String,Any}("J1_meV" => 0.25, "J2_meV" => 0.01, "gzz" => 3.35,
                      "sigma_gzz" => 0.8, "sigma_J" => 0.5)
params, _ = SV.sv_apply_param_overrides(p0, Dict("param_overrides" => ov))

cuts = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
cut = cuts[3]        # (0.5,0,0) at 9 T
ctx = SV.sv_kpm_context(params, controls; field_T = cut.field_T, realization = 0)

# Warm the KPM path so the timed loop is compute, not JIT.
warm = SV.sv_kpm_spectrum_from_context(ctx, controls, cut)
say("$(TAG) warmup: nq=$(warm.q_samples) chunks=$(warm.n_chunks) ",
    @sprintf("%.1f s  threads=%d", warm.seconds, Threads.nthreads()))

# Rendezvous so all processes compute simultaneously.
bdir = getenv("MP_BARRIER", "")
if !isempty(bdir)
    mkpath(bdir); touch(joinpath(bdir, "ready_$(TAG)"))
    while count(startswith("ready_"), readdir(bdir)) < NPROC
        sleep(0.05)
    end
end

t_start = time()
secs = Float64[]
for i in 1:NSPEC
    s = SV.sv_kpm_spectrum_from_context(ctx, controls, cut)
    push!(secs, s.seconds)
end
t_end = time()

rss_mb = try
    parse(Int, split(read("/proc/self/status", String))[findfirst(==("VmRSS:"),
          split(read("/proc/self/status", String))) + 1]) / 1024
catch; NaN end

say(@sprintf("RESULT\t%s\tnspec=%d\tnq=%d\tchunks=%d\tcompute_s=%.2f\tspan_s=%.2f\tper_spec_s=%.2f\tt0=%.3f\tt1=%.3f\trss_MB=%.0f",
             TAG, NSPEC, warm.q_samples, warm.n_chunks, sum(secs), t_end - t_start,
             sum(secs) / NSPEC, t_start, t_end, rss_mb))
say("MP_DONE $(TAG)")
