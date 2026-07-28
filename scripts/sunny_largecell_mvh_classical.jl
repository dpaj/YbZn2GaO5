#!/usr/bin/env julia

# Classical finite-temperature M(H,T) for the MINIMAL single-disordered-phase
# Sunny model, in the same effective-Hamiltonian scheme as the KPM LSWT neutron
# calculation (same P1 triangular net, same J1/J2 shells, same sigma_J bond
# disorder, same per-site sigma_gzz, same field direction, same seed).
#
# What is deliberately different from the analytical co-fit:
#
#   * ONE disordered phase. No separate non-dispersive/"flat" component. The
#     premise under test is that thermal population of the low-energy modes that
#     disorder generates inside a single coupled Hamiltonian does the job that
#     the analytical model does with a second independent phase.
#
#   * Finite temperature by classical Boltzmann sampling rather than T = 0
#     minimization.
#
# Known systematic, quantified by the validation panel: classical statistics
# assigns every magnon mode an occupancy kT/eps where quantum statistics gives
# exp(-eps/kT). Above saturation the quantum answer is M_sat = gzz*S to ~1e-8,
# so the Langevin deficit measured there is pure classical over-counting. On this
# Hamiltonian at 0.42 K it is about 0.11 uB/site (~6%) at 9 T, and it is
# field-dependent, so it distorts the SHAPE of M(H) and not merely the scale.
#
# Run with:
#   julia --project=. scripts/sunny_largecell_mvh_classical.jl

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using Random
using CairoMakie
using Sunny

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

const SV = SunnyValidation
const KB_MEV_PER_K = 0.08617333262   # matches SV_KB_MEV_PER_K

# -----------------------------------------------------------------------------
# Config plumbing (same pattern as scripts/sunny_kpm_1d_disp_grid_2sigmaJ_vs_exp.jl)
# -----------------------------------------------------------------------------

function _deepmerge!(dest::Dict, src::Dict)
    for (k, v) in src
        if v isa Dict && haskey(dest, k) && dest[k] isa Dict
            _deepmerge!(dest[k], v)
        else
            dest[k] = v
        end
    end
    return dest
end

function _repo_path(repo_root::AbstractString, p::AbstractString)
    isabspath(p) && return normpath(p)
    return normpath(joinpath(repo_root, splitpath(p)...))
end

const DEFAULT_DIAG_CONTROLS = "configs/sunny_largecell_mvh_classical_controls.toml"

# Config path resolution order: first positional ARG, then SUNNY_MVH_CONTROLS,
# then the default. Lets a cheap smoke config be swapped in without editing the
# production one.
function _diag_controls_path(repo_root::AbstractString)
    if !isempty(ARGS) && !isempty(ARGS[1])
        return _repo_path(repo_root, ARGS[1])
    end
    env = get(ENV, "SUNNY_MVH_CONTROLS", "")
    isempty(env) || return _repo_path(repo_root, env)
    return _repo_path(repo_root, DEFAULT_DIAG_CONTROLS)
end

function _load_diagnostic_controls(repo_root::AbstractString)
    diag_path = _diag_controls_path(repo_root)
    isfile(diag_path) || error("Could not find diagnostic controls: $diag_path")
    diag = load_toml_config(diag_path)

    base_rel = get(get(diag, "paths", Dict{String,Any}()), "base_controls_toml",
                   "configs/sunny_validation_controls.toml")
    base_path = _repo_path(repo_root, String(base_rel))
    controls = load_toml_config(base_path)

    haskey(diag, "control_overrides") && _deepmerge!(controls, diag["control_overrides"])

    if haskey(diag, "paths")
        paths = diag["paths"]
        haskey(paths, "figure_subdir") && (controls["paths"]["figure_subdir"] = paths["figure_subdir"])
        haskey(paths, "table_subdir") && (controls["paths"]["table_subdir"] = paths["table_subdir"])
    end

    return (; diag, controls, diag_path, base_path)
end

_tuple3(v) = (Int(v[1]), Int(v[2]), Int(v[3]))

function _repeat_factor(cell_size, seed_dims)
    cs, sd = _tuple3(cell_size), _tuple3(seed_dims)
    rf = ntuple(i -> cs[i] ÷ sd[i], 3)
    for i in 1:3
        rf[i] * sd[i] == cs[i] ||
            error("cell_size $cs is not an integer multiple of seed_dims $sd along axis $i")
    end
    return rf
end

# -----------------------------------------------------------------------------
# Timing / profile bookkeeping
# -----------------------------------------------------------------------------

function _timed!(rows::Vector{NamedTuple}, stage, f; cell="", realization=-1, sampler="", field_T=NaN, note="")
    local val
    t = @elapsed begin
        val = f()
    end
    push!(rows, (; stage=String(stage), cell=String(cell), realization=Int(realization),
        sampler=String(sampler), field_T=Float64(field_T), seconds=Float64(t), note=String(note)))
    return val
end

function _write_profile_csv(path, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "stage,cell,realization,sampler,field_T,seconds,note")
        for r in rows
            println(io, join(Any[r.stage, r.cell, r.realization, r.sampler,
                @sprintf("%.6g", r.field_T), @sprintf("%.9g", r.seconds),
                replace(r.note, "," => ";")], ","))
        end
    end
    return path
end

# -----------------------------------------------------------------------------
# System construction: minimal single disordered phase
# -----------------------------------------------------------------------------

function _build_system(params, controls::Dict, cell_size, seed_dims, realization::Int;
                       field_T::Real=0.0)
    seed_dims_t = _tuple3(seed_dims)
    rf = _repeat_factor(cell_size, seed_dims_t)
    lc = get(controls, "largecell", Dict{String,Any}())

    base = SV.sv_build_effective_sunny_system(params, controls;
        component=:dispersive, dims=seed_dims_t, field_T=field_T)
    sys = rf == (1, 1, 1) ? to_inhomogeneous(base.sys) :
                            to_inhomogeneous(repeat_periodically(base.sys, rf))

    SV.sv_apply_disorder!(sys, params, controls;
        component=:dispersive,
        include_exchange=Bool(get(lc, "include_exchange_disorder", true)),
        include_gzz=Bool(get(lc, "include_gzz_disorder", true)),
        realization=realization)

    return (; sys, crystal=base.crystal, units=base.units, repeat_factor=rf)
end

function _initialize_field_polarized!(sys, uhat; field_T::Real=1.0)
    dip = collect((Float64(field_T) < 0 ? -1.0 : 1.0) .* uhat)
    for site in eachsite(sys)
        set_dipole!(sys, dip, site)
    end
    return sys
end

_set_field!(sys, uhat, units, B_T) = set_field!(sys, collect(uhat .* (Float64(B_T) * units.T)))

# -----------------------------------------------------------------------------
# Observable. sv_m_parallel_uB_per_site returns mean(g*S . uhat); Sunny's moment
# is -g*S, so moment_sign (-1 by convention in [largecell]) flips it into the
# experimental positive-M convention.
# -----------------------------------------------------------------------------

_M_site(sys, uhat, moment_sign) = moment_sign * SV.sv_m_parallel_uB_per_site(sys, uhat)

# Fully saturated moment along uhat for THIS disorder realization. Per-site gzz
# is disordered, so mean(g_i)*S differs from gzz*S by O(sigma_gzz/sqrt(N)) and a
# small cell can legitimately exceed the nominal gzz*S. Deficits must therefore
# be measured against this, not against gzz*S.
function _M_sat_site(sys, uhat, spin_S)
    acc = 0.0
    n = 0
    for site in eachsite(sys)
        acc += dot(uhat, sys.gs[site] * uhat) * spin_S
        n += 1
    end
    return n > 0 ? acc / n : NaN
end

function _langevin_dt(diag_run::Dict, sys, integrator_kT, damping)
    mode = String(get(diag_run, "dt_mode", "manual"))
    dt_manual = Float64(get(diag_run, "dt", 0.02))
    mode == "manual" && return (dt_manual, "manual")
    tol = Float64(get(diag_run, "suggest_timestep_tol", 0.01))
    safety = Float64(get(diag_run, "dt_safety", 0.25))
    try
        probe = Langevin(dt_manual; damping=damping, kT=integrator_kT)
        bound = Sunny.suggest_timestep_aux(sys, probe; tol=tol)
        return (safety * Float64(bound), @sprintf("suggest bound=%.6g safety=%.3g", bound, safety))
    catch err
        @warn "suggest_timestep_aux unavailable; falling back to manual dt" error=err dt=dt_manual
        return (dt_manual, "manual_fallback")
    end
end

# Returns (M_mean, M_std_of_samples).
function _sample_M(sys, integrator, uhat, moment_sign; n_sample::Int, stride::Int)
    vals = Vector{Float64}(undef, n_sample)
    for i in 1:n_sample
        for _ in 1:stride
            step!(sys, integrator)
        end
        vals[i] = _M_site(sys, uhat, moment_sign)
    end
    return (mean(vals), n_sample > 1 ? std(vals) : 0.0)
end

function _thermalize!(sys, integrator, nsteps::Int)
    for _ in 1:nsteps
        step!(sys, integrator)
    end
    return sys
end

# -----------------------------------------------------------------------------
# One field sweep for a given (cell size, realization, sampler)
# -----------------------------------------------------------------------------

function _sweep(params, controls::Dict, diag_run::Dict, cell_size, seed_dims, realization::Int,
                sampler::AbstractString, Bs, profile::Vector{NamedTuple})
    uhat = SV.sv_field_direction(controls)
    units = SV.sv_units()
    moment_sign = Float64(get(get(controls, "largecell", Dict{String,Any}()), "moment_sign", -1.0))
    T_K = Float64(get(diag_run, "temperature_K", 0.42))
    kT = KB_MEV_PER_K * T_K
    maxiters = Int(get(diag_run, "minimize_maxiters", 10_000))
    damping = Float64(get(diag_run, "damping", 0.1))
    n_th = Int(get(diag_run, "n_thermalize", 3000))
    n_sa = Int(get(diag_run, "n_sample", 500))
    stride = Int(get(diag_run, "sample_stride", 10))
    rescale = Bool(get(diag_run, "spin_rescaling_for_static_sum_rule", false))
    adiabatic = Bool(get(diag_run, "adiabatic_field_continuation", true))
    relax_before_thermalize = Bool(get(diag_run, "relax_before_thermalize", true))
    init_state = String(get(diag_run, "initial_spin_state", "field_polarized"))
    cell_label = join(string.(_tuple3(cell_size)), "x")

    built = _timed!(profile, "build_system", () ->
        _build_system(params, controls, cell_size, seed_dims, realization; field_T=0.0);
        cell=cell_label, realization=realization, sampler=sampler,
        note=@sprintf("seed_dims=%s", join(string.(_tuple3(seed_dims)), "x")))
    sys = built.sys
    nsites = prod(_tuple3(cell_size))

    if rescale
        set_spin_rescaling_for_static_sum_rule!(sys)
    end

    if init_state == "field_polarized"
        _initialize_field_polarized!(sys, uhat; field_T=1.0)
    elseif init_state == "random"
        Random.seed!(Int(controls["common"]["seed"]) + 101 + 7919 * realization)
        randomize_spins!(sys)
    end

    dt, dt_note = _langevin_dt(diag_run, sys, kT, damping)

    M = zeros(Float64, length(Bs))
    Msd = zeros(Float64, length(Bs))

    for (i, B) in enumerate(Bs)
        _set_field!(sys, uhat, units, B)
        if !adiabatic
            _initialize_field_polarized!(sys, uhat; field_T=1.0)
        end

        if sampler == "minimize_energy"
            _timed!(profile, "minimize_energy", () -> minimize_energy!(sys; maxiters=maxiters);
                cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                note=@sprintf("maxiters=%d", maxiters))
            M[i] = _M_site(sys, uhat, moment_sign)
            Msd[i] = 0.0

        elseif sampler == "langevin"
            # Relax to the classical ground state at this field FIRST, then heat
            # to kT. Thermalizing straight from the previous field's thermal
            # state does not converge in a practical number of steps: the smoke
            # run without this step returned M(7 T) ~ 0.02 uB instead of ~1.73.
            if relax_before_thermalize
                _timed!(profile, "minimize_energy", () -> minimize_energy!(sys; maxiters=maxiters);
                    cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                    note=@sprintf("pre-thermalize maxiters=%d", maxiters))
            end
            integ = Langevin(dt; damping=damping, kT=kT)
            _timed!(profile, "thermalize", () -> _thermalize!(sys, integ, n_th);
                cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                note=@sprintf("dt=%.6g damping=%.3g kT=%.6g n=%d [%s]", dt, damping, kT, n_th, dt_note))
            m, s = _timed!(profile, "sample", () ->
                _sample_M(sys, integ, uhat, moment_sign; n_sample=n_sa, stride=stride);
                cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                note=@sprintf("n_sample=%d stride=%d", n_sa, stride))
            M[i], Msd[i] = m, s

        elseif sampler == "langevin_then_midpoint"
            if relax_before_thermalize
                _timed!(profile, "minimize_energy", () -> minimize_energy!(sys; maxiters=maxiters);
                    cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                    note=@sprintf("pre-thermalize maxiters=%d", maxiters))
            end
            integ = Langevin(dt; damping=damping, kT=kT)
            _timed!(profile, "thermalize", () -> _thermalize!(sys, integ, n_th);
                cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                note=@sprintf("dt=%.6g damping=%.3g kT=%.6g n=%d [%s]", dt, damping, kT, n_th, dt_note))
            imp = ImplicitMidpoint(dt)
            m, s = _timed!(profile, "sample_microcanonical", () ->
                _sample_M(sys, imp, uhat, moment_sign; n_sample=n_sa, stride=stride);
                cell=cell_label, realization=realization, sampler=sampler, field_T=B,
                note=@sprintf("n_sample=%d stride=%d", n_sa, stride))
            M[i], Msd[i] = m, s

        else
            error("Unknown sampler '$sampler'. Use minimize_energy, langevin, or langevin_then_midpoint.")
        end
    end

    M_sat_real = _M_sat_site(sys, uhat, Float64(controls["common"]["spin_S"]))
    return (; M, Msd, cell_label, nsites, repeat_factor=built.repeat_factor, dt, dt_note, kT,
              M_sat_realization=M_sat_real)
end

# -----------------------------------------------------------------------------
# Validation above saturation: quantum answer is M_sat = gzz*S to ~1e-8 there,
# so the Langevin deficit is pure classical over-counting.
# -----------------------------------------------------------------------------

function _lswt_magnon_depletion(sys, kT; min_gap_over_kT::Real=2.0)
    # Field-polarized collinear state only: M_z = g*(S - <n_i>), so the
    # per-site depletion is sum_nu n_B(eps_nu) / N. This is MEANINGLESS below
    # saturation, where the spectrum is gapless and sum n_B diverges — it can
    # return a "depletion" larger than the total moment. `valid` reports whether
    # the gap is comfortably above kT; callers must check it.
    #
    # Strong gzz disorder pushes saturation well past the nominal clean B_sat,
    # because low-g sites need a much larger field, so the validation fields
    # have to be chosen against the DISORDERED saturation field, not the clean one.
    swt = SpinWaveTheory(sys; measure=nothing)
    epsv = vec(dispersion(swt, [[0.0, 0.0, 0.0]]))
    epsv = epsv[isfinite.(epsv)]
    nb = 0.0
    for e in epsv
        e > 1e-9 || continue
        x = e / kT
        x < 700 && (nb += 1.0 / expm1(x))
    end
    gap = isempty(epsv) ? NaN : minimum(epsv)
    valid = isfinite(gap) && gap > min_gap_over_kT * kT
    return (; n_total=nb, n_per_site=nb / max(1, length(epsv)),
              eps_min=gap, eps_max=isempty(epsv) ? NaN : maximum(epsv),
              nmodes=length(epsv), valid=valid)
end

function _run_validation(params, controls::Dict, diag_run::Dict, cell_size, seed_dims,
                         profile::Vector{NamedTuple})
    uhat = SV.sv_field_direction(controls)
    units = SV.sv_units()
    moment_sign = Float64(get(get(controls, "largecell", Dict{String,Any}()), "moment_sign", -1.0))
    T_K = Float64(get(diag_run, "temperature_K", 0.42))
    kT = KB_MEV_PER_K * T_K
    fields = Float64.(get(diag_run, "validation_fields_T", [9.0, 14.0]))
    do_lswt = Bool(get(diag_run, "validation_lswt", true))
    maxiters = Int(get(diag_run, "minimize_maxiters", 10_000))
    damping = Float64(get(diag_run, "damping", 0.1))
    n_th = Int(get(diag_run, "n_thermalize", 3000))
    n_sa = Int(get(diag_run, "n_sample", 500))
    stride = Int(get(diag_run, "sample_stride", 10))
    cell_label = join(string.(_tuple3(cell_size)), "x")
    spin_S = Float64(controls["common"]["spin_S"])

    rows = NamedTuple[]
    for B in fields
        built = _build_system(params, controls, cell_size, seed_dims, 0; field_T=B)
        sys = built.sys
        # Re-assert the field explicitly rather than trusting it to survive
        # repeat_periodically / to_inhomogeneous.
        _set_field!(sys, uhat, units, B)
        _initialize_field_polarized!(sys, uhat; field_T=B)
        minimize_energy!(sys; maxiters=maxiters)
        M_T0 = _M_site(sys, uhat, moment_sign)
        M_sat = _M_sat_site(sys, uhat, spin_S)

        lswt = do_lswt ? _timed!(profile, "lswt_validation", () -> _lswt_magnon_depletion(sys, kT);
                    cell=cell_label, sampler="lswt", field_T=B) :
                    (; n_total=NaN, n_per_site=NaN, eps_min=NaN, eps_max=NaN, nmodes=-1, valid=false)
        # Only meaningful when the magnon spectrum is gapped; otherwise the Bose
        # sum diverges and the number would be garbage rather than small.
        dM_lswt = lswt.valid ? params.gzz * lswt.n_per_site : NaN

        dt, _ = _langevin_dt(diag_run, sys, kT, damping)
        integ = Langevin(dt; damping=damping, kT=kT)
        _thermalize!(sys, integ, n_th)
        M_lang, M_lang_sd = _sample_M(sys, integ, uhat, moment_sign; n_sample=n_sa, stride=stride)

        # M_sat here is mean(g_i)*S for this realization, not the nominal gzz*S.
        # Above saturation the quantum answer is M_sat to ~1e-8 (see dM_lswt_bose),
        # so M_sat - M_langevin is the classical over-counting systematic.
        overcount = M_sat - M_lang
        saturated = abs(M_sat - M_T0) < 1e-3 * max(eps(), abs(M_sat))
        push!(rows, (; field_T=B, cell=cell_label, M_sat_realization=M_sat, M_T0=M_T0,
            M_langevin=M_lang, M_langevin_std=M_lang_sd,
            T0_deficit_vs_sat_uB=M_sat - M_T0,
            eps_min_meV=lswt.eps_min, eps_max_meV=lswt.eps_max, nmodes=lswt.nmodes,
            gap_over_kT=lswt.eps_min / kT, lswt_valid=lswt.valid, T0_saturated=saturated,
            dM_lswt_bose=dM_lswt,
            classical_overcount_uB=overcount,
            classical_overcount_percent=100 * overcount / max(eps(), abs(M_sat))))
        @printf("  validation B=%5.2f T : M_sat(real)=%.6f  M_T0=%.6f  M_Langevin=%.6f  classical over-count=%.6f uB (%.2f%%)\n",
                B, M_sat, M_T0, M_lang, overcount,
                100 * overcount / max(eps(), abs(M_sat)))
        if lswt.valid && saturated
            @printf("                       gap=%.4g meV (%.1f kT)  LSWT+Bose dM=%.3g uB  -> over-count is the whole deficit\n",
                    lswt.eps_min, lswt.eps_min / kT, dM_lswt)
        else
            @printf("                       NOT a valid calibration point: gap=%.4g meV (%.2f kT), T=0 state %s saturated. Raise validation_fields_T.\n",
                    lswt.eps_min, lswt.eps_min / kT, saturated ? "is" : "is NOT")
        end
    end
    return rows
end

# -----------------------------------------------------------------------------
# dt / damping convergence scan. Equilibrium averages must be independent of
# damping; if they are not, dt is too large.
# -----------------------------------------------------------------------------

function _run_convergence_scan(params, controls::Dict, diag_run::Dict, cell_size, seed_dims)
    dts = Float64.(get(diag_run, "dt_scan", Float64[]))
    dampings = Float64.(get(diag_run, "damping_scan", Float64[]))
    (isempty(dts) || isempty(dampings)) && return NamedTuple[]

    uhat = SV.sv_field_direction(controls)
    units = SV.sv_units()
    moment_sign = Float64(get(get(controls, "largecell", Dict{String,Any}()), "moment_sign", -1.0))
    T_K = Float64(get(diag_run, "temperature_K", 0.42))
    kT = KB_MEV_PER_K * T_K
    B = Float64(get(diag_run, "convergence_field_T", 9.0))
    maxiters = Int(get(diag_run, "minimize_maxiters", 10_000))
    n_th = Int(get(diag_run, "n_thermalize", 3000))
    n_sa = Int(get(diag_run, "n_sample", 500))
    stride = Int(get(diag_run, "sample_stride", 10))
    cell_label = join(string.(_tuple3(cell_size)), "x")

    rows = NamedTuple[]
    println("  dt / damping convergence at B = $(B) T (equilibrium M must not depend on damping):")
    for dt in dts, damping in dampings
        built = _build_system(params, controls, cell_size, seed_dims, 0; field_T=B)
        sys = built.sys
        _set_field!(sys, uhat, units, B)
        _initialize_field_polarized!(sys, uhat; field_T=B)
        minimize_energy!(sys; maxiters=maxiters)
        integ = Langevin(dt; damping=damping, kT=kT)
        _thermalize!(sys, integ, n_th)
        m, s = _sample_M(sys, integ, uhat, moment_sign; n_sample=n_sa, stride=stride)
        push!(rows, (; cell=cell_label, field_T=B, dt=dt, damping=damping,
            M=m, M_std=s, kT=kT))
        @printf("    dt=%-6.4g damping=%-5.3g : M = %+.6f +- %.6f uB\n", dt, damping, m, s)
    end
    return rows
end

# -----------------------------------------------------------------------------
# Scale handling
# -----------------------------------------------------------------------------

function _best_positive_scale(y_exp, x_model)
    mask = isfinite.(y_exp) .& isfinite.(x_model)
    any(mask) || return 1.0
    num = sum(y_exp[mask] .* x_model[mask])
    den = sum(x_model[mask] .^ 2)
    den > 0 || return 1.0
    return max(0.0, num / den)
end

# Joint nonnegative least squares of y ~ a*x1 + b*x2.
#
# Used to fit the moment amplitude A_M and the Van Vleck slope TOGETHER, because
# chi_vv_muB_per_T was fitted inside the analytical model and carries that
# model's assumptions. Holding it fixed here can manufacture a residual that
# looks like missing physics. Here x1 = the Sunny moment and x2 = B, so
# A_M = a and the effective Van Vleck slope is b/a.
function _best_two_component_scale(y, x1, x2)
    m = isfinite.(y) .& isfinite.(x1) .& isfinite.(x2)
    any(m) || return (1.0, 0.0)
    yv, u, v = y[m], x1[m], x2[m]
    a11, a12, a22 = sum(u .^ 2), sum(u .* v), sum(v .^ 2)
    b1, b2 = sum(yv .* u), sum(yv .* v)
    det = a11 * a22 - a12^2
    if abs(det) < 1e-300
        return (_best_positive_scale(yv, u), 0.0)
    end
    a = (b1 * a22 - b2 * a12) / det
    b = (b2 * a11 - b1 * a12) / det
    if a >= 0 && b >= 0
        return (a, b)
    end
    # Unconstrained optimum is outside the nonnegative quadrant: compare the two
    # single-component boundary solutions and keep the better residual.
    aa = _best_positive_scale(yv, u)
    bb = _best_positive_scale(yv, v)
    return sum((yv .- aa .* u) .^ 2) <= sum((yv .- bb .* v) .^ 2) ? (aa, 0.0) : (0.0, bb)
end

# Overrides let a by-eye or exploratory parameter set be tested without touching
# configs/best_fit_parameters.toml. Every override is echoed and recorded.
function _apply_param_overrides(params, run::Dict)
    ov = get(run, "param_overrides", Dict{String,Any}())
    (ov isa Dict && !isempty(ov)) || return (params, String[])
    applied = String[]
    kv = Dict{Symbol,Float64}()
    for k in sort(collect(keys(ov)))
        sym = Symbol(k)
        hasproperty(params, sym) ||
            error("[run.param_overrides] key '$k' is not a canonical model parameter")
        push!(applied, @sprintf("%s: %.6g -> %.6g",
            k, Float64(getproperty(params, sym)), Float64(ov[k])))
        kv[sym] = Float64(ov[k])
    end
    return (merge(params, NamedTuple(kv)), applied)
end

# "clean" = identical Hamiltonian with both disorder widths set to zero, so
# sv_apply_disorder! writes uniform J and uniform gzz.
_variant_params(params, variant::AbstractString) =
    variant == "clean" ? merge(params, (; sigma_J=0.0, sigma_gzz=0.0)) : params

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

function _write_results_csv(path, curves, Bs, exp_interp, diag_run, A_M_fixed, free_scales, joint_scales)
    mkpath(dirname(path))
    include_vv = Bool(get(diag_run, "include_chi_vv", true))
    open(path, "w") do io
        println(io, join([
            "cell", "nsites", "variant", "realization", "sampler", "B_T",
            "M_disp_raw_uB_per_site", "M_disp_std_uB_per_site",
            "M_vv_uB", "M_model_unscaled_uB", "M_total_fixed_scale_uB", "M_total_free_scale_uB",
            "M_total_joint_fit_uB", "M_vv_joint_fit_uB",
            "M_exp_interp_uB_per_Yb", "residual_fixed_scale_uB", "residual_joint_fit_uB",
            "magnetization_global_scale_fixed", "magnetization_global_scale_free",
            "magnetization_global_scale_joint", "chi_vv_joint_fit_muB_per_T",
            "temperature_K", "kT_meV", "dt", "damping", "n_thermalize", "n_sample", "sample_stride",
            "include_chi_vv", "chi_vv_muB_per_T", "gzz", "J1_meV", "J2_meV",
            "sigma_J", "sigma_gzz", "spin_S", "seed", "M_sat_realization_uB", "model",
        ], ","))
        for c in curves
            key = (c.cell_label, c.sampler, c.variant)
            A_free = get(free_scales, key, NaN)
            aj, bj = get(joint_scales, key, (NaN, NaN))
            chi_j = aj > 0 ? bj / aj : NaN
            for (i, B) in enumerate(Bs)
                vv = include_vv ? c.params_used.chi_vv_muB_per_T * B : 0.0
                unscaled = c.M[i] + vv
                tot_fixed = A_M_fixed * unscaled
                tot_free = A_free * unscaled
                tot_joint = aj * c.M[i] + bj * B
                println(io, join(Any[
                    c.cell_label, c.nsites, c.variant, c.realization, c.sampler, @sprintf("%.10g", B),
                    @sprintf("%.10g", c.M[i]), @sprintf("%.10g", c.Msd[i]),
                    @sprintf("%.10g", vv), @sprintf("%.10g", unscaled),
                    @sprintf("%.10g", tot_fixed), @sprintf("%.10g", tot_free),
                    @sprintf("%.10g", tot_joint), @sprintf("%.10g", bj * B),
                    @sprintf("%.10g", exp_interp[i]), @sprintf("%.10g", tot_fixed - exp_interp[i]),
                    @sprintf("%.10g", tot_joint - exp_interp[i]),
                    @sprintf("%.10g", A_M_fixed), @sprintf("%.10g", A_free),
                    @sprintf("%.10g", aj), @sprintf("%.10g", chi_j),
                    @sprintf("%.10g", c.T_K), @sprintf("%.10g", c.kT),
                    @sprintf("%.10g", c.dt), @sprintf("%.10g", c.damping),
                    c.n_thermalize, c.n_sample, c.sample_stride,
                    include_vv, @sprintf("%.10g", c.params_used.chi_vv_muB_per_T),
                    @sprintf("%.10g", c.params_used.gzz), @sprintf("%.10g", c.params_used.J1_meV),
                    @sprintf("%.10g", c.params_used.J2_meV), @sprintf("%.10g", c.params_used.sigma_J),
                    @sprintf("%.10g", c.params_used.sigma_gzz), c.spin_S, c.seed,
                    @sprintf("%.10g", c.M_sat_realization),
                    "minimal_single_disordered_phase",
                ], ","))
            end
        end
    end
    return path
end

function _write_manifest_csv(path, curves, A_M_fixed, free_scales, joint_scales, diag_run,
                            sunny_version, overrides_applied)
    mkpath(dirname(path))
    ov = isempty(overrides_applied) ? "none" : join(overrides_applied, "; ")
    open(path, "w") do io
        println(io, join([
            "cell", "nsites", "variant", "seed_dims", "repeat_factor", "realization", "sampler",
            "nB", "Bmin_T", "Bmax_T", "temperature_K", "kT_meV",
            "dt", "dt_note", "damping", "n_thermalize", "n_sample", "sample_stride",
            "magnetization_global_scale_fixed", "magnetization_global_scale_free",
            "magnetization_global_scale_joint", "chi_vv_joint_fit_muB_per_T",
            "M_at_max_field_fixed_scale", "M_sat_realization_uB",
            "gzz", "J1_meV", "J2_meV", "sigma_J", "sigma_gzz", "chi_vv_muB_per_T",
            "param_overrides", "spin_rescaling_for_static_sum_rule",
            "include_chi_vv", "model", "sunny_version",
        ], ","))
        for c in curves
            key = (c.cell_label, c.sampler, c.variant)
            A_free = get(free_scales, key, NaN)
            aj, bj = get(joint_scales, key, (NaN, NaN))
            include_vv = Bool(get(diag_run, "include_chi_vv", true))
            vv_end = include_vv ? c.params_used.chi_vv_muB_per_T * c.Bmax : 0.0
            println(io, join(Any[
                c.cell_label, c.nsites, c.variant, c.seed_dims_label, c.repeat_label,
                c.realization, c.sampler,
                c.nB, @sprintf("%.10g", c.Bmin), @sprintf("%.10g", c.Bmax),
                @sprintf("%.10g", c.T_K), @sprintf("%.10g", c.kT),
                @sprintf("%.10g", c.dt), replace(c.dt_note, "," => ";"),
                @sprintf("%.10g", c.damping), c.n_thermalize, c.n_sample, c.sample_stride,
                @sprintf("%.10g", A_M_fixed), @sprintf("%.10g", A_free),
                @sprintf("%.10g", aj), @sprintf("%.10g", aj > 0 ? bj / aj : NaN),
                @sprintf("%.10g", A_M_fixed * (c.M[end] + vv_end)),
                @sprintf("%.10g", c.M_sat_realization),
                @sprintf("%.10g", c.params_used.gzz), @sprintf("%.10g", c.params_used.J1_meV),
                @sprintf("%.10g", c.params_used.J2_meV), @sprintf("%.10g", c.params_used.sigma_J),
                @sprintf("%.10g", c.params_used.sigma_gzz),
                @sprintf("%.10g", c.params_used.chi_vv_muB_per_T),
                replace(ov, "," => ";"),
                Bool(get(diag_run, "spin_rescaling_for_static_sum_rule", false)),
                include_vv, "minimal_single_disordered_phase", string(sunny_version),
            ], ","))
        end
    end
    return path
end

function _write_simple_csv(path, header::Vector{String}, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join([v isa AbstractString ? replace(v, "," => ";") :
                              (v isa Integer ? string(v) : @sprintf("%.10g", Float64(v)))
                              for v in values(r)], ","))
        end
    end
    return path
end

function _make_plot(path, curves, Bs, data, A_M_fixed, joint_scales, diag_run)
    mkpath(dirname(path))
    ylim = get(diag_run, "plot_ylim", nothing)
    suffix = String(get(diag_run, "title_suffix", "minimal single disordered phase"))
    T_K = Float64(get(diag_run, "temperature_K", 0.42))

    cells = unique([c.cell_label for c in curves])
    sort!(cells; by=cl -> first(c.nsites for c in curves if c.cell_label == cl))
    prod_cell = cells[end]
    variants = unique([c.variant for c in curves])
    samp = "langevin" in [c.sampler for c in curves] ? "langevin" :
           first(unique([c.sampler for c in curves]))
    vcolor = Dict("disordered" => :darkorange, "clean" => :royalblue)

    _mean_raw(sel) = vec(mean(reduce(hcat, [c.M for c in sel]); dims=2))

    fig = Figure(size=(1350, 560))

    # Panel 1: disordered versus clean at the joint (A_M, chi_vv) fit.
    ax1 = Axis(fig[1, 1]; xlabel="B (T)", ylabel="M (μB / Yb)", xticks=0:1:ceil(Bs[end]),
        title=@sprintf("%s cell, %s at %.2f K — A_M and χ_vv fitted jointly", prod_cell, samp, T_K))
    scatter!(ax1, data.B_T, data.M_muB_per_Yb; markersize=3, color=:black, label="experiment")
    for v in variants
        sel = [c for c in curves if c.variant == v && c.sampler == samp && c.cell_label == prod_cell]
        isempty(sel) && continue
        aj, bj = get(joint_scales, (prod_cell, samp, v), (A_M_fixed, 0.0))
        raw = _mean_raw(sel)
        col = get(vcolor, v, :seagreen)
        lines!(ax1, Bs, aj .* raw .+ bj .* Bs; color=col, linewidth=2,
            label=@sprintf("%s (A_M=%.3f, χ_vv=%.3f)", v, aj, aj > 0 ? bj / aj : NaN))
        lines!(ax1, Bs, bj .* Bs; color=col, linestyle=:dot, linewidth=1.5,
            label=@sprintf("%s Van Vleck part", v))
    end
    ylim !== nothing && length(ylim) == 2 && ylims!(ax1, Float64(ylim[1]), Float64(ylim[2]))
    axislegend(ax1; position=:lt, framevisible=false, labelsize=9)

    # Panel 2: residual of the joint fit — this is where a missing component shows.
    ax2 = Axis(fig[1, 2]; xlabel="B (T)", ylabel="model − experiment (μB / Yb)",
        xticks=0:1:ceil(Bs[end]),
        title="Residual at joint fit (all cells, all variants)")
    hlines!(ax2, [0.0]; color=:black, linestyle=:dash, linewidth=1)
    exp_on_grid = Float64.(SV.sv_interp1(data.B_T, data.M_muB_per_Yb, Bs))
    for v in variants, (k, cl) in enumerate(cells)
        sel = [c for c in curves if c.variant == v && c.sampler == samp && c.cell_label == cl]
        isempty(sel) && continue
        aj, bj = get(joint_scales, (cl, samp, v), (A_M_fixed, 0.0))
        res = aj .* _mean_raw(sel) .+ bj .* Bs .- exp_on_grid
        keep = isfinite.(res)
        lines!(ax2, Bs[keep], res[keep]; color=get(vcolor, v, :seagreen),
            linestyle=(:solid, :dash, :dot)[min(k, 3)], linewidth=2,
            label=@sprintf("%s %s (max |r|=%.3f)", v, cl, maximum(abs, res[keep])))
    end
    axislegend(ax2; position=:lb, framevisible=false, labelsize=9)

    Label(fig[0, :], "Sunny large-cell classical M(H,T): " * suffix; fontsize=15, font=:bold)
    save(path, fig)
    return path
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

function main()
    (; diag, controls, diag_path, base_path) = _load_diagnostic_controls(REPO_ROOT)
    run = get(diag, "run", Dict{String,Any}())

    (; params, path) = SV.sv_load_params(REPO_ROOT, controls)
    print_canonical_model_parameters(params)

    # A_M is read BEFORE overrides so the fixed reference stays the canonical
    # analytical value even when the Hamiltonian parameters are overridden.
    A_M_fixed = SV.sv_magnetization_global_scale(params, controls)

    params, overrides_applied = _apply_param_overrides(params, run)
    if !isempty(overrides_applied)
        println()
        println("[run.param_overrides] applied — these are NOT the canonical best-fit values:")
        for s in overrides_applied
            println("  ", s)
        end
    end

    T_K = Float64(get(run, "temperature_K", 0.42))
    kT = KB_MEV_PER_K * T_K
    Bmin = Float64(get(run, "Bmin_T", 0.0))
    Bmax = Float64(get(run, "Bmax_T", 7.0))
    nB = Int(get(run, "nB", 36))
    Bs = collect(range(Bmin, Bmax; length=nB))
    seed_dims = get(run, "seed_dims", [3, 3, 1])
    cell_sizes = get(run, "cell_sizes", [[36, 36, 1]])
    n_real = Int(get(run, "n_realizations", 3))
    samplers = String.(get(run, "samplers", ["minimize_energy", "langevin"]))
    variants = String.(get(run, "disorder_variants", ["disordered"]))
    M_sat = params.gzz * Float64(controls["common"]["spin_S"])

    @info "Sunny large-cell classical M(H,T): minimal single disordered phase" sunny_version=SV.sv_try_pkgversion(Sunny) params_path=path diag_path base_path
    @printf("\nT = %.4g K -> kT = %.6g meV ;  J1 = %.6g meV ;  M_sat = gzz*S = %.6f uB/site\n",
            T_K, kT, params.J1_meV, M_sat)
    @printf("Field grid: %.3g to %.3g T, %d points.  Cell sizes: %s.  Realizations: %d.\n",
            Bmin, Bmax, nB, join([join(string.(_tuple3(c)), "x") for c in cell_sizes], ", "), n_real)
    @printf("Samplers: %s.  Disorder variants: %s.\n", join(samplers, ", "), join(variants, ", "))
    @printf("sigma_J = %.6g (fractional), sigma_gzz = %.6g, chi_vv = %.6g uB/T\n",
            params.sigma_J, params.sigma_gzz, params.chi_vv_muB_per_T)
    @printf("Fixed magnetization_global_scale A_M = %.10g\n\n", A_M_fixed)

    mag_path = SV.sv_repo_path(REPO_ROOT, controls["paths"]["magnetization_csv"])
    data = SV.sv_read_magnetization_csv(mag_path)
    exp_interp = Float64.(SV.sv_interp1(data.B_T, data.M_muB_per_Yb, Bs))

    profile = NamedTuple[]
    curves = NamedTuple[]

    for cell_size in cell_sizes
        cell_label = join(string.(_tuple3(cell_size)), "x")
        for variant in variants
            pv = _variant_params(params, variant)
            # A clean system has no disorder, so every realization is identical.
            nr = variant == "clean" ? 1 : n_real
            for sampler in samplers
                for r in 0:(nr - 1)
                    @printf("--- cell=%s variant=%-11s sampler=%-22s realization=%d\n",
                            cell_label, variant, sampler, r)
                    sw = _sweep(pv, controls, run, cell_size, seed_dims, r, sampler, Bs, profile)
                    push!(curves, (;
                        cell_label=sw.cell_label, nsites=sw.nsites, realization=r, sampler=sampler,
                        variant=variant, params_used=pv,
                        M=sw.M, Msd=sw.Msd, dt=sw.dt, dt_note=sw.dt_note, kT=sw.kT, T_K=T_K,
                        M_sat_realization=sw.M_sat_realization,
                    damping=Float64(get(run, "damping", 0.1)),
                    n_thermalize=Int(get(run, "n_thermalize", 3000)),
                    n_sample=Int(get(run, "n_sample", 500)),
                    sample_stride=Int(get(run, "sample_stride", 10)),
                    nB=nB, Bmin=Bmin, Bmax=Bmax,
                    seed_dims_label=join(string.(_tuple3(seed_dims)), "x"),
                    repeat_label=join(string.(sw.repeat_factor), "x"),
                        spin_S=Float64(controls["common"]["spin_S"]),
                        seed=Int(controls["common"]["seed"]),
                    ))
                    @printf("    M(%.3g T) = %+.6f uB/site (raw), fixed-scale total = %+.6f uB/Yb\n",
                            Bmax, sw.M[end],
                            A_M_fixed * (sw.M[end] + (Bool(get(run, "include_chi_vv", true)) ?
                                pv.chi_vv_muB_per_T * Bmax : 0.0)))
                end
            end
        end
    end

    # Scales per (cell, sampler, variant), averaged over disorder realizations.
    #   free  : amplitude only, chi_vv held at the analytical value
    #   joint : amplitude AND Van Vleck slope fitted together
    include_vv = Bool(get(run, "include_chi_vv", true))
    win = Float64.(get(run, "scale_fit_window_T", [Bmin, Bmax]))
    inwin = (Bs .>= win[1]) .& (Bs .<= win[2])
    keys3 = unique([(c.cell_label, c.sampler, c.variant) for c in curves])
    free_scales = Dict{Tuple{String,String,String},Float64}()
    joint_scales = Dict{Tuple{String,String,String},Tuple{Float64,Float64}}()
    for k in keys3
        sel = [c for c in curves
               if c.cell_label == k[1] && c.sampler == k[2] && c.variant == k[3]]
        isempty(sel) && continue
        raw = vec(mean(reduce(hcat, [c.M for c in sel]); dims=2))
        vv = include_vv ? sel[1].params_used.chi_vv_muB_per_T .* Bs : zeros(length(Bs))
        free_scales[k] = _best_positive_scale(exp_interp[inwin], (raw .+ vv)[inwin])
        joint_scales[k] = _best_two_component_scale(exp_interp[inwin], raw[inwin], Bs[inwin])
    end

    println()
    @printf("Scale fits (fixed reference A_M = %.6f, analytical chi_vv = %.6f uB/T):\n",
            A_M_fixed, params.chi_vv_muB_per_T)
    @printf("  %-10s %-11s %-22s %10s %10s %12s\n",
            "cell", "variant", "sampler", "A_M_free", "A_M_joint", "chi_vv_joint")
    for k in sort(keys3)
        aj, bj = joint_scales[k]
        @printf("  %-10s %-11s %-22s %10.5f %10.5f %12.5f\n",
                k[1], k[3], k[2], free_scales[k], aj, aj > 0 ? bj / aj : NaN)
    end
    println()

    validation_rows = NamedTuple[]
    if Bool(get(run, "validate_above_saturation", true))
        # Calibrate the classical over-counting on the CLEAN system by default.
        # Strong gzz disorder pushes the disordered saturation field far above the
        # experimental range (low-g sites need many times more field), so a
        # disordered cell is not gapped at 9-14 T and the LSWT reference is invalid
        # there. The per-mode classical over-counting is a generic property of
        # classical statistics, so the clean system calibrates it cleanly.
        vv_name = String(get(run, "validation_variant", "clean"))
        println("Validation above saturation on the '$vv_name' system (quantum M = M_sat when gapped; deficit is classical over-counting):")
        validation_rows = _run_validation(_variant_params(params, vv_name), controls, run,
                                         cell_sizes[end], seed_dims, profile)
        println()
    end

    conv_rows = _run_convergence_scan(params, controls, run, cell_sizes[end], seed_dims)
    isempty(conv_rows) || println()

    out_tab = SV.sv_repo_path(REPO_ROOT, controls["paths"]["table_subdir"])
    out_fig = SV.sv_repo_path(REPO_ROOT, controls["paths"]["figure_subdir"])
    mkpath(out_tab); mkpath(out_fig)

    csv_path = joinpath(out_tab, "sunny_largecell_mvh_classical.csv")
    man_path = joinpath(out_tab, "sunny_largecell_mvh_classical_manifest.csv")
    prof_path = joinpath(out_tab, "sunny_largecell_mvh_classical_profile.csv")
    fig_path = joinpath(out_fig, "sunny_largecell_mvh_classical.png")

    _write_results_csv(csv_path, curves, Bs, exp_interp, run, A_M_fixed, free_scales, joint_scales)
    _write_manifest_csv(man_path, curves, A_M_fixed, free_scales, joint_scales, run,
                        SV.sv_try_pkgversion(Sunny), overrides_applied)
    _write_profile_csv(prof_path, profile)
    _make_plot(fig_path, curves, Bs, data, A_M_fixed, joint_scales, run)

    val_path = ""
    if !isempty(validation_rows)
        val_path = joinpath(out_tab, "sunny_largecell_mvh_classical_validation.csv")
        _write_simple_csv(val_path, String.(collect(keys(validation_rows[1]))), validation_rows)
    end
    conv_path = ""
    if !isempty(conv_rows)
        conv_path = joinpath(out_tab, "sunny_largecell_mvh_classical_convergence.csv")
        _write_simple_csv(conv_path, String.(collect(keys(conv_rows[1]))), conv_rows)
    end

    println("Wrote:")
    for p in [csv_path, man_path, prof_path, val_path, conv_path, fig_path]
        isempty(p) || println("  ", p)
    end
    @printf("\nTotal wall clock in profiled stages: %.1f s\n", sum(r -> r.seconds, profile; init=0.0))

    return (; Bs, curves, free_scales, joint_scales, A_M_fixed, overrides_applied,
              validation_rows, conv_rows,
              csv_path, man_path, prof_path, fig_path, val_path, conv_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
