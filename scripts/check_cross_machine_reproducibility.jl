#!/usr/bin/env julia
# Cross-machine reproducibility of the neutron objective, against a target measured on the DGX.
#
#   julia -t auto --project=. scripts/check_cross_machine_reproducibility.jl
#
# TARGET: sv_neutron_objective chi2_red = 27.047850, six cuts, 81 q, n = 8.
#
# Measured on neutrons-dgx01 at repo 2717e10 with Julia 1.12.6, Sunny 0.9.1
# (git-tree-sha1 12f5d5415d334b2d0f28a9be39771902d70c7847), 32 threads, BLAS pinned to 1.
# That box also reproduced it identically to six decimals across 8/16/32/64 threads, across
# 4-way concurrent chains, and across runs days apart -- so it is bit-reproducible WITHIN a
# machine. This checks the harder claim: ACROSS machines, with a different OS, a different core
# count, a different chunk count and a different Julia patch version.
#
# WHY IT IS WORTH ONE EVALUATION. Every parameter comparison the two machines have exchanged --
# chi2 surfaces, gzz determinations, per-cut breakdowns -- assumes the objective means the same
# thing on both. If it does not, those comparisons are silently incommensurable, and that is the
# kind of error that is worthless to find later and cheap to find now. Nelder-Mead also depends
# on a deterministic objective, so this is a precondition for the fit rather than a nicety.
#
# Every control is set EXPLICITLY here rather than inherited, because several differ from the
# base config -- notably n_energy = 241 and the energy range 0.1-4.2 meV, where the base config
# has 161. Inheriting would silently test a different configuration.
#
# On a mismatch, compare E/site per (field, realization) FIRST: that isolates the ground state
# from the KPM path, and a different local minimum under a different BLAS/LAPACK build is the
# most likely cause. The script prints it for exactly that reason.

using Printf, Statistics, LinearAlgebra, Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl")); using .YZGOCofit
include(joinpath(REPO_ROOT, "src", "sunny_validation.jl")); using .SunnyValidation
const SV = SunnyValidation

const TARGET = 27.047850
BLAS.set_num_threads(1)          # as on the DGX; BLAS threads inside a threaded q loop contend

controls = SV.sv_load_controls(REPO_ROOT)

kc = controls["kpm"]
kc["dims"] = [3, 3, 1]
kc["system_size"] = [36, 36, 1]
kc["method"] = "lanczos"
kc["tol"] = 0.05
kc["maxiters"] = 1000
kc["regularization"] = 1e-5
kc["energy_min_meV"] = 0.1
kc["energy_max_meV"] = 4.2
kc["n_energy"] = 241
kc["kernel_fwhm_meV"] = 0.05
kc["thread_max_chunks"] = 32
kc["thread_min_q_per_chunk"] = 1
kc["thread_min_q"] = 1
kc["qtags"] = ["0_1_0", "0p33_0p33_0", "0p5_0_0"]

eh = kc["experimental_histogram"]
eh["enabled"] = true
eh["mode"] = "analytical_cut_volume_grid"
eh["measured_grid_mode"] = "uniform_grid"
eh["n_measured_h"] = 3; eh["n_measured_k"] = 3; eh["n_measured_l"] = 1

qa = kc["q_averaging"]
qa["enabled"] = true
qa["mode"] = "gaussian_grid"
qa["n_h"] = 3; qa["n_k"] = 3; qa["n_l"] = 1
qa["sigma_H_rlu"] = 0.037; qa["sigma_K_rlu"] = 0.037; qa["sigma_L_rlu"] = 0.0
# resolution_quadrature deliberately left at its default, gauss_hermite (d2166ff). Setting it
# would defeat the purpose of checking that the DEFAULT agrees across machines.

cm = controls["common"]
cm["seed"] = 20260611
cm["spin_S"] = 0.5
cm["field_direction"] = [0.0, 0.0, 1.0]

(; params) = SV.sv_load_params(REPO_ROOT, controls)
params = merge(params, (; J1_meV=0.15, J2_meV=0.01, gzz=3.50, sigma_gzz=0.80, sigma_J=0.50))

cuts = SV.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
nq = length(SV.sv_kpm_1d_q_sampler(cuts[1], controls).qs)

@printf("environment: julia %s, %d threads, BLAS %d, Sunny %s\n",
        VERSION, Threads.nthreads(), BLAS.get_num_threads(),
        string(SV.sv_try_pkgversion(Sunny)))
@printf("chunks for nq = %d: %d   (DGX used 32; a different value is part of the test)\n",
        nq, SV.sv_kpm_q_chunks(controls, nq))
@printf("cuts: %d  %s\n", length(cuts),
        join(["$(c.qtag)@$(Int(c.field_T))T" for c in cuts], " "))
@printf("params: J1=%.4f J2=%.4f gzz=%.4f sigma_gzz=%.4f sigma_J=%.4f\n\n",
        params.J1_meV, params.J2_meV, params.gzz, params.sigma_gzz, params.sigma_J)

t = time()
o = SV.sv_neutron_objective(params, controls, cuts; realizations=0:7, threaded=true,
        maxiters=1000, on_failure=:record)
el = time() - t

# Ground state first: it isolates a relaxation difference from a KPM difference.
println("E/site per (field, realization) -- compare these FIRST on a mismatch:")
for c in sort(o.contexts; by = x -> (x.field_T, x.realization))
    @printf("  B=%5.2f T  r=%d  E/site = %.10f  converged=%s  maxiters=%d\n",
            c.field_T, c.realization, c.E_per_site, c.converged, c.maxiters)
end

println("\nper-cut chi2_red:")
for pc in o.per_cut
    @printf("  %-14s %5.1f T  %.6f\n", pc.qtag, pc.field_T, pc.chi2_red)
end

d = o.chi2_red - TARGET
@printf("\nchi2_red = %.6f   target = %.6f   difference = %+.3e  (%.2e relative)\n",
        o.chi2_red, TARGET, d, abs(d) / TARGET)
@printf("scale = %.8g   regularization used = %s   failures = %d   (%.0f s, KPM %.0f s)\n",
        o.scale, string(o.regularization_values), o.n_failed, el, o.kpm_seconds)

if abs(d) < 5e-6
    println("\nMATCH to six decimals. The objective is platform-independent across OS, core " *
            "count, chunk count and Julia patch version, so every cross-machine number is " *
            "comparable and the optimizer has the determinism it needs.")
elseif abs(d) / TARGET < 1e-3
    println("\nCLOSE but not bit-identical. Something differs at the round-off level -- most " *
            "likely a BLAS/LAPACK build difference in minimize_energy!. Compare E/site above " *
            "against the DGX's before looking anywhere else.")
else
    println("\nMISMATCH. Do NOT treat cross-machine chi2 comparisons as meaningful until this " *
            "is resolved. Check, in order: relax_attempts (must be 1), the quadrature default, " *
            "then E/site per context to separate the ground state from the KPM path.")
end
