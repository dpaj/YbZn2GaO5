#!/usr/bin/env julia

# Fixed-model Sunny KPM GPU summary for YbZn2GaO5.
#
# This is a science-facing diagnostic script, not an optimizer.
# It loads the existing YZGO repo controls/data, evaluates the fixed best-fit
# dispersive Sunny KPM model for the configured 1D neutron cuts, and makes a
# combined neutron + mean-field magnetization summary figure.
#
# Recommended GPU-branch use from repo root:
#
#   $env:SUNNY_KPM_ENABLE_GPU="1"
#   $env:SUNNY_KPM_GPU_BACKEND="CUDA"
#   $env:SUNNY_KPM_GPU_BATCHED="1"
#   $env:SUNNY_KPM_GPU_PRECISION="Float64"
#   julia --project=envs\sunny-kpm-gpu scripts\sunny_kpm_gpu_fixed_model_summary.jl
#
# CPU fallback:
#
#   $env:SUNNY_KPM_ENABLE_GPU="0"
#   julia --project=. scripts\sunny_kpm_gpu_fixed_model_summary.jl
#
# Warm interactive use:
#
#   julia --project=envs\sunny-kpm-gpu
#
#   include("scripts/sunny_kpm_gpu_fixed_model_summary.jl")
#   run_sunny_kpm_gpu_fixed_model_summary(; enable_gpu=true,
#                                         backend="CUDA",
#                                         batched=true,
#                                         precision="Float64",
#                                         Ei_meV=[4.65],
#                                         neutron_intensity_ylim=(0.0, 0.003),
#                                         include_flat_component=false,
#                                         kpm_energy_max_meV=4.2,
#                                         kpm_n_energy=241,
#                                         energy_resample_mode="direct_interpolation",
#                                         param_overrides=(; J1_meV=0.25,
#                                                           J2_meV=0.01,
#                                                           sigma_J=0.30,
#                                                           gzz=3.8,
#                                                           sigma_gzz=0.4))
#
# Re-run the same function in the same Julia session for warm timings. The script
# intentionally does not auto-run when included interactively.

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation

using Sunny
using CairoMakie
using LinearAlgebra
using Printf
using Statistics
using Dates
using TOML

# -----------------------------------------------------------------------------
# Small config / IO helpers
# -----------------------------------------------------------------------------

const CONTROL_PATH = joinpath(REPO_ROOT, "configs", "sunny_kpm_gpu_fixed_model_controls.toml")
const INTERACTIVE_WORKFLOW_OVERRIDES = Ref{Dict{String,Any}}(Dict{String,Any}())
const INTERACTIVE_PARAM_OVERRIDES = Ref{Any}(nothing)

function _as_bool(x, default=false)
    x === nothing && return default
    x isa Bool && return x
    s = lowercase(string(x))
    return s in ("1", "true", "yes", "on")
end

function _recursive_merge!(dst::Dict, src::Dict)
    for (k, v) in src
        if haskey(dst, k) && dst[k] isa Dict && v isa Dict
            _recursive_merge!(dst[k], v)
        else
            dst[k] = v
        end
    end
    return dst
end

function _apply_interactive_control_overrides!(controls::Dict, diag_controls::Dict)
    wf = get(diag_controls, "workflow", Dict{String,Any}())
    kpm_overrides = Dict{String,Any}()

    haskey(wf, "kpm_energy_min_meV") && (kpm_overrides["energy_min_meV"] = Float64(wf["kpm_energy_min_meV"]))
    haskey(wf, "kpm_energy_max_meV") && (kpm_overrides["energy_max_meV"] = Float64(wf["kpm_energy_max_meV"]))
    haskey(wf, "kpm_n_energy") && (kpm_overrides["n_energy"] = Int(wf["kpm_n_energy"]))
    haskey(wf, "kpm_kernel_fwhm_meV") && (kpm_overrides["kernel_fwhm_meV"] = Float64(wf["kpm_kernel_fwhm_meV"]))
    haskey(wf, "energy_resample_mode") && (kpm_overrides["energy_resample_mode"] = String(wf["energy_resample_mode"]))

    if !isempty(kpm_overrides)
        get!(controls, "kpm", Dict{String,Any}())
        _recursive_merge!(controls["kpm"], kpm_overrides)
    end
    return controls
end

function _load_workflow_controls(repo_root::AbstractString)
    base_controls = SunnyValidation.sv_load_controls(repo_root)
    diag_controls = isfile(CONTROL_PATH) ? TOML.parsefile(CONTROL_PATH) : Dict{String,Any}()

    if haskey(diag_controls, "control_overrides")
        _recursive_merge!(base_controls, diag_controls["control_overrides"])
    end

    # Environment variables are convenient for the DGX/A100 GPU workflow.
    gpu = get!(diag_controls, "gpu", Dict{String,Any}())
    gpu["enabled"] = _as_bool(get(ENV, "SUNNY_KPM_ENABLE_GPU", get(gpu, "enabled", false)))
    gpu["backend"] = get(ENV, "SUNNY_KPM_GPU_BACKEND", get(gpu, "backend", "CUDA"))
    gpu["batched"] = _as_bool(get(ENV, "SUNNY_KPM_GPU_BATCHED", get(gpu, "batched", true)))
    gpu["precision"] = get(ENV, "SUNNY_KPM_GPU_PRECISION", get(gpu, "precision", "Float64"))

    # Interactive wrapper overrides. These intentionally affect only the
    # diagnostic workflow layer, not the canonical repo controls.
    overrides = INTERACTIVE_WORKFLOW_OVERRIDES[]
    if !isempty(overrides)
        wf = get!(diag_controls, "workflow", Dict{String,Any}())
        for (k, v) in overrides
            wf[k] = v
        end
    end

    _apply_interactive_control_overrides!(base_controls, diag_controls)

    return base_controls, diag_controls
end

function _workflow_array(diag::Dict, key::AbstractString, fallback)
    wf = get(diag, "workflow", Dict{String,Any}())
    return haskey(wf, key) ? collect(wf[key]) : collect(fallback)
end

function _workflow_get(diag::Dict, key::AbstractString, fallback)
    wf = get(diag, "workflow", Dict{String,Any}())
    return get(wf, key, fallback)
end

function _workflow_optional_pair(diag::Dict, key::AbstractString)
    wf = get(diag, "workflow", Dict{String,Any}())
    haskey(wf, key) || return nothing
    v = collect(wf[key])
    length(v) == 2 || error("workflow.$key must have two values, e.g. [0.0, 0.003]")
    return (Float64(v[1]), Float64(v[2]))
end

function _output_dirs(repo_root::AbstractString, diag::Dict)
    subdir = String(_workflow_get(diag, "output_subdir", "sunny_validation/kpm_gpu_fixed_model"))
    table_dir = joinpath(repo_root, "results", "feature_tables", splitpath(subdir)...)
    fig_dir = joinpath(repo_root, "results", "figures", splitpath(subdir)...)
    mkpath(table_dir)
    mkpath(fig_dir)
    return table_dir, fig_dir
end

function _namedtuple_from_dict(d::AbstractDict)
    names = Tuple(Symbol.(collect(keys(d))))
    vals = Tuple(collect(values(d)))
    return NamedTuple{names}(vals)
end

function _coerce_param_overrides(x)
    x === nothing && return nothing
    x isa NamedTuple && return x
    x isa AbstractDict && return _namedtuple_from_dict(x)
    if x isa AbstractVector{<:Pair}
        return _namedtuple_from_dict(Dict(Symbol(k) => v for (k, v) in x))
    end
    error("param_overrides must be a NamedTuple, Dict, or vector of Pair.")
end

function _apply_param_overrides(params)
    overrides = INTERACTIVE_PARAM_OVERRIDES[]
    overrides === nothing && return params
    override_nt = _coerce_param_overrides(overrides)
    return merge(params, override_nt)
end

function _csv_cell(x)
    if x isa Real
        return isfinite(Float64(x)) ? @sprintf("%.10g", Float64(x)) : string(x)
    elseif x === missing || x === nothing
        return ""
    else
        s = string(x)
        if occursin(",", s) || occursin("\"", s) || occursin("\n", s)
            return "\"" * replace(s, "\"" => "\"\"") * "\""
        else
            return s
        end
    end
end

function _write_rows_csv(path::AbstractString, rows::Vector{<:NamedTuple})
    mkpath(dirname(path))
    isempty(rows) && error("No rows to write for $path")
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for r in rows
            println(io, join((_csv_cell(getproperty(r, nm)) for nm in names), ","))
        end
    end
    return path
end

# -----------------------------------------------------------------------------
# Lightweight timers
# -----------------------------------------------------------------------------

mutable struct TimerTable
    rows::Vector{NamedTuple}
end
TimerTable() = TimerTable(NamedTuple[])

function _timed!(f, timers::TimerTable; stage::String, Ei_meV=NaN, field_T=NaN,
                 qtag="", component="", q_samples=0, note="")
    t0 = time_ns()
    val = f()
    seconds = (time_ns() - t0) / 1e9
    push!(timers.rows, (; stage, Ei_meV=Float64(Ei_meV), field_T=Float64(field_T),
                         qtag=String(qtag), component=String(component),
                         q_samples=Int(q_samples), seconds=seconds,
                         seconds_per_q=q_samples > 0 ? seconds / q_samples : NaN,
                         note=String(note)))
    @printf("%-22s %8.3f s", stage, seconds)
    q_samples > 0 && @printf("   (%d Q, %.6f s/Q)", q_samples, seconds / q_samples)
    isempty(note) || @printf("   %s", note)
    println()
    return val
end

# -----------------------------------------------------------------------------
# Optional GPU helpers
# -----------------------------------------------------------------------------

function _require_pkg(uuid::String, name::String)
    return Base.require(Base.PkgId(Base.UUID(uuid), name))
end

function _gpu_precision_type(diag::Dict)
    s = lowercase(String(get(get(diag, "gpu", Dict{String,Any}()), "precision", "Float64")))
    s in ("float32", "f32", "single") && return Float32
    s in ("float64", "f64", "double") && return Float64
    error("Unsupported GPU precision=$s. Use Float32 or Float64.")
end

function _maybe_make_backend(diag::Dict)
    gpu = get(diag, "gpu", Dict{String,Any}())
    _as_bool(get(gpu, "enabled", false)) || return nothing

    backend_name = lowercase(String(get(gpu, "backend", "CUDA")))
    if backend_name == "cuda"
        try
            _require_pkg("63c18a36-062a-441e-b654-da1e3ab1ce7c", "KernelAbstractions")
            cuda = _require_pkg("052768ef-5323-5732-b1bb-66c8b64840ba", "CUDA")
            functional = Base.invokelatest(getproperty(cuda, :functional))
            if !functional
                @warn "CUDA loaded but CUDA.functional() is false; using CPU fallback"
                return nothing
            end
            CUDABackend = getproperty(cuda, :CUDABackend)
            return Base.invokelatest(CUDABackend)
        catch err
            @warn "Could not initialize CUDA backend; using CPU fallback" exception=(err, catch_backtrace())
            return nothing
        end
    else
        @warn "Only CUDA is wired in this first fixed-model GPU workflow; using CPU fallback" backend=backend_name
        return nothing
    end
end

function _make_device_kpm(swt, backend, diag::Dict)
    gpu = get(diag, "gpu", Dict{String,Any}())
    batched = _as_bool(get(gpu, "batched", true))
    precision = _gpu_precision_type(diag)
    if batched && isdefined(Sunny, :to_device_batched)
        return Base.invokelatest(getproperty(Sunny, :to_device_batched), swt, backend; precision=precision)
    elseif isdefined(Sunny, :to_device)
        return Base.invokelatest(getproperty(Sunny, :to_device), swt, backend)
    else
        error("Sunny GPU extension entry points not found. Expected Sunny.to_device_batched or Sunny.to_device.")
    end
end


# -----------------------------------------------------------------------------
# Energy-axis resampling helpers
# -----------------------------------------------------------------------------

function _interp1_zero(x::AbstractVector, y::AbstractVector, xout::AbstractVector)
    n = length(x)
    n == length(y) || error("_interp1_zero: x and y lengths differ")
    n >= 2 || error("_interp1_zero: need at least two source grid points")
    out = zeros(Float64, length(xout))
    x1 = Float64(first(x))
    xN = Float64(last(x))
    for (i, xo0) in pairs(xout)
        xo = Float64(xo0)
        if xo < x1 || xo > xN
            out[i] = 0.0
        else
            j = searchsortedlast(x, xo)
            if j <= 0
                out[i] = Float64(y[1])
            elseif j >= n
                out[i] = Float64(y[end])
            else
                xa = Float64(x[j])
                xb = Float64(x[j+1])
                t = xb == xa ? 0.0 : (xo - xa) / (xb - xa)
                out[i] = (1 - t) * Float64(y[j]) + t * Float64(y[j+1])
            end
        end
    end
    return out
end

function _model_to_experimental_energy_grid(model_energy, model_intensity, exp_energy, controls::Dict)
    kc = controls["kpm"]
    mode = Symbol(get(kc, "energy_resample_mode", "repo_resolved"))

    if mode in (:repo_resolved, :sv_resolved, :sunny_validation)
        hist_mode = Symbol(get(kc, "histogram_mode", "bin_average"))
        return SunnyValidation.sv_model_to_experimental_energy_grid_resolved(
            model_energy, model_intensity, exp_energy, controls; section="kpm", mode=hist_mode
        )
    elseif mode in (:direct_interpolation, :interp, :linear_interp)
        return _interp1_zero(model_energy, model_intensity, exp_energy)
    else
        error("Unsupported kpm.energy_resample_mode=$mode. Use repo_resolved or direct_interpolation.")
    end
end

# -----------------------------------------------------------------------------
# KPM context and spectra
# -----------------------------------------------------------------------------

function _initialize_field_polarized!(sys, controls::Dict)
    u = SunnyValidation.sv_field_direction(controls)
    for site in eachsite(sys)
        set_dipole!(sys, u, site)
    end
    return sys
end

function _construct_kpm(sys, controls::Dict)
    kc = controls["kpm"]
    measure = SunnyValidation.sv_sunny_measure(sys, controls)
    tol = Float64(kc["tol"])
    method = Symbol(get(kc, "method", "lanczos"))
    try
        return SpinWaveTheoryKPM(sys; measure=measure, tol=tol, method=method)
    catch err
        @warn "SpinWaveTheoryKPM(...; method=$method) failed; retrying without method keyword" exception=(err, catch_backtrace())
        return SpinWaveTheoryKPM(sys; measure=measure, tol=tol)
    end
end

function _prepare_kpm_context(params, controls::Dict, diag::Dict, timers::TimerTable;
                              component::Symbol, field_T::Real, backend=nothing)
    kc = controls["kpm"]
    sizectl = SunnyValidation.sv_system_size_controls(controls, "kpm")
    dims = sizectl.dims
    repeat_factor = sizectl.repeat_factor
    include_exchange = get(kc, "include_exchange_disorder", true)
    include_gzz = get(kc, "include_gzz_disorder", true)
    maxiters = Int(kc["maxiters"])

    sys = _timed!(timers; stage="build_system", field_T, component=String(component)) do
        base = SunnyValidation.sv_build_effective_sunny_system(params, controls; component, dims, field_T)
        s = base.sys
        if repeat_factor != (1, 1, 1)
            s = to_inhomogeneous(repeat_periodically(s, repeat_factor))
        else
            s = to_inhomogeneous(s)
        end
        SunnyValidation.sv_apply_disorder!(s, params, controls; component,
                                           include_exchange=include_exchange,
                                           include_gzz=include_gzz)
        return s
    end

    init_mode = Symbol(_workflow_get(diag, "initial_spin_state", "field_polarized"))
    _timed!(timers; stage="initialize_spins", field_T, component=String(component)) do
        if init_mode == :field_polarized
            _initialize_field_polarized!(sys, controls)
        elseif init_mode == :random
            randomize_spins!(sys)
        else
            error("Unsupported initial_spin_state=$init_mode. Use field_polarized or random.")
        end
    end

    if _as_bool(_workflow_get(diag, "relax_ground_state", true))
        _timed!(timers; stage="minimize_energy", field_T, component=String(component)) do
            minimize_energy!(sys; maxiters=maxiters)
        end
    end

    swt = _timed!(timers; stage="construct_kpm", field_T, component=String(component)) do
        _construct_kpm(sys, controls)
    end

    eval_obj = swt
    gpu_used = false
    if backend !== nothing
        eval_obj = _timed!(timers; stage="transfer_to_device", field_T,
                           component=String(component),
                           note="GPU transfer") do
            _make_device_kpm(swt, backend, diag)
        end
        gpu_used = true
    end

    kc = controls["kpm"]
    energies = collect(range(Float64(kc["energy_min_meV"]),
                             Float64(kc["energy_max_meV"]);
                             length=Int(kc["n_energy"])))
    kernel = gaussian(fwhm=Float64(kc["kernel_fwhm_meV"]))

    return (; sys, swt, eval_obj, energies, kernel, gpu_used, component,
            field_T=Float64(field_T), system_size=sizectl.system_size)
end

function _intensities(eval_obj, qs, energies, kernel)
    # Use invokelatest because the KA extension methods may load after script compilation.
    try
        return Base.invokelatest(Sunny.intensities, eval_obj, qs;
                                 energies=energies, kernel=kernel, kT=0.0, verbose=false)
    catch err
        # Stable CPU Sunny may not accept verbose/kT for this method combination.
        return Base.invokelatest(Sunny.intensities, eval_obj, qs;
                                 energies=energies, kernel=kernel)
    end
end

function _component_cut_spectrum_from_context(ctx, params, controls::Dict,
                                              cut, timers::TimerTable)
    sampler = SunnyValidation.sv_kpm_1d_q_sampler(cut, controls)
    qs = sampler.qs
    res = _timed!(timers; stage="kpm_intensities",
                  Ei_meV=cut.Ei_meV, field_T=cut.field_T,
                  qtag=cut.qtag, component=String(ctx.component),
                  q_samples=length(qs),
                  note=ctx.gpu_used ? "GPU" : "CPU") do
        _intensities(ctx.eval_obj, qs, ctx.energies, ctx.kernel)
    end

    raw = SunnyValidation.sv_try_extract_sunny_intensity(res)
    I0 = SunnyValidation.sv_orient_sunny_intensity_matrix(raw, length(ctx.energies), length(qs))
    Ipost, ff_w, ff_a, qmag = SunnyValidation.sv_apply_form_factor_to_intensity(I0, qs, controls)
    Iavg = SunnyValidation.sv_kpm_1d_average_qsampled_intensity(Ipost, sampler)

    Imodel_expgrid = _model_to_experimental_energy_grid(ctx.energies, Iavg, cut.energy_meV, controls)

    return (; energy_model_meV=ctx.energies,
            intensity_model_grid=Iavg,
            energy_exp_meV=cut.energy_meV,
            intensity_exp=cut.intensity,
            intensity_model_expgrid=Imodel_expgrid,
            q_samples=length(qs),
            q_measured_samples=sampler.n_measured,
            q_resolution_samples=sampler.n_resolution,
            form_factor_weight=sum(sampler.weights .* ff_w),
            form_factor_amplitude=sum(sampler.weights .* ff_a),
            qmag_Ainv=sum(sampler.weights .* qmag),
            gpu_used=ctx.gpu_used)
end

function _panel_scale(model::AbstractVector, data::AbstractVector, params, controls::Dict, diag::Dict)
    mode = Symbol(_workflow_get(diag, "neutron_scale_mode", "panel_least_squares"))
    if mode == :panel_least_squares
        return SunnyValidation.sv_best_vertical_scale(model, data)
    elseif mode == :best_fit
        return SunnyValidation.sv_neutron_scale(params, controls)
    elseif mode == :none || mode == :unit
        return 1.0
    else
        error("Unsupported workflow.neutron_scale_mode=$mode")
    end
end

function _controls_for_cut(base_controls::Dict; Ei_meV::Real, field_T::Real, qtag::AbstractString)
    c = deepcopy(base_controls)
    c["kpm"]["Ei_meV"] = Float64(Ei_meV)
    c["kpm"]["qtags"] = [String(qtag)]
    c["common"]["fields_T"] = [Float64(field_T)]
    return c
end

function _try_load_cut(repo_root, base_controls; Ei_meV, field_T, qtag)
    c = _controls_for_cut(base_controls; Ei_meV, field_T, qtag)
    try
        cuts = SunnyValidation.sv_load_kpm_experimental_cuts(repo_root, c)
        if isempty(cuts)
            return nothing, c
        end
        return first(cuts), c
    catch err
        @warn "Skipping missing/problematic neutron cut" Ei_meV field_T qtag exception=(err, catch_backtrace())
        return nothing, c
    end
end

# -----------------------------------------------------------------------------
# Plotting and summary CSV
# -----------------------------------------------------------------------------

function _plot_neutron_panel!(ax, item, diag::Dict)
    cut = item.cut
    lines!(ax, item.energy_meV, item.model_scaled, linewidth=2, label="Sunny KPM")
    scatter!(ax, item.energy_meV, item.intensity_exp, markersize=6, label="experiment")
    ax.title = @sprintf("Ei %.2f meV, %.0f T, %s", cut.Ei_meV, cut.field_T,
                        SunnyValidation.sv_qtag_label(cut.qtag))
    ax.xlabel = "Energy transfer (meV)"
    ax.ylabel = "Intensity (arb.)"
    ylim = _workflow_optional_pair(diag, "neutron_intensity_ylim")
    ylim === nothing || ylims!(ax, ylim[1], ylim[2])
    return ax
end

function _write_neutron_long_csv(path, items)
    rows = NamedTuple[]
    for item in items
        cut = item.cut
        for i in eachindex(item.energy_meV)
            push!(rows, (;
                Ei_meV=cut.Ei_meV,
                field_T=cut.field_T,
                qtag=cut.qtag,
                energy_meV=item.energy_meV[i],
                intensity_exp=item.intensity_exp[i],
                intensity_model_unscaled=item.model_unscaled[i],
                intensity_model_scaled=item.model_scaled[i],
                panel_scale=item.panel_scale,
                q_samples=item.q_samples,
                q_measured_samples=item.q_measured_samples,
                q_resolution_samples=item.q_resolution_samples,
                gpu_used=item.gpu_used ? 1 : 0,
                component_mode=item.component_mode,
                energy_resample_mode=item.energy_resample_mode,
                form_factor_weight=item.form_factor_weight,
                form_factor_amplitude=item.form_factor_amplitude,
                qmag_Ainv=item.qmag_Ainv,
            ))
        end
    end
    isempty(rows) || _write_rows_csv(path, rows)
    return path
end

function _write_neutron_raw_energy_csv(path, items)
    rows = NamedTuple[]
    for item in items
        cut = item.cut
        for i in eachindex(item.energy_model_meV)
            push!(rows, (;
                Ei_meV=cut.Ei_meV,
                field_T=cut.field_T,
                qtag=cut.qtag,
                energy_model_meV=item.energy_model_meV[i],
                intensity_model_raw_unscaled=item.model_raw_unscaled[i],
                intensity_model_raw_scaled=item.model_raw_scaled[i],
                panel_scale=item.panel_scale,
                q_samples=item.q_samples,
                q_measured_samples=item.q_measured_samples,
                q_resolution_samples=item.q_resolution_samples,
                gpu_used=item.gpu_used ? 1 : 0,
                component_mode=item.component_mode,
                energy_resample_mode=item.energy_resample_mode,
            ))
        end
    end
    isempty(rows) || _write_rows_csv(path, rows)
    return path
end

function _plot_summary(fig_path, neutron_items, mag_result, mag_data, Eis, fields, qtags, diag::Dict)
    nrows_neutron = length(Eis) * length(fields)
    ncols = length(qtags)
    fig = Figure(size=(360 * ncols, 260 * nrows_neutron + 360))

    item_by_key = Dict{Tuple{Float64,Float64,String},Any}()
    for item in neutron_items
        cut = item.cut
        item_by_key[(cut.Ei_meV, cut.field_T, cut.qtag)] = item
    end

    row = 1
    for Ei in Eis
        for B in fields
            for (j, qtag) in enumerate(qtags)
                ax = Axis(fig[row, j])
                key = (Float64(Ei), Float64(B), String(qtag))
                if haskey(item_by_key, key)
                    _plot_neutron_panel!(ax, item_by_key[key], diag)
                    if row == 1 && j == ncols
                        axislegend(ax, position=:rt)
                    end
                else
                    ax.title = @sprintf("Ei %.2f meV, %.0f T, %s", Float64(Ei), Float64(B),
                                        SunnyValidation.sv_qtag_label(qtag))
                    text!(ax, 0.5, 0.5, text="no cut loaded", align=(:center, :center), space=:relative)
                    ylim = _workflow_optional_pair(diag, "neutron_intensity_ylim")
                    ylim === nothing || ylims!(ax, ylim[1], ylim[2])
                    hidedecorations!(ax)
                end
            end
            row += 1
        end
    end

    axm = Axis(fig[row, 1:ncols], xlabel="B (T)", ylabel="M (μB / Yb)",
               title="Mean-field bridge magnetization")
    scatter!(axm, mag_data.B_T, mag_data.M_muB_per_Yb, markersize=6, label="experiment")
    lines!(axm, mag_result.B_T, mag_result.M_total, linewidth=2, label="total")
    lines!(axm, mag_result.B_T, mag_result.M_disp, linestyle=:dash, label="dispersive")
    lines!(axm, mag_result.B_T, mag_result.M_flat, linestyle=:dashdot, label="flat")
    lines!(axm, mag_result.B_T, mag_result.M_vv, linestyle=:dot, label="Van Vleck")
    axislegend(axm, position=:rb)

    save(fig_path, fig)
    return fig_path
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

function main()
    controls, diag = _load_workflow_controls(REPO_ROOT)
    params_loaded = SunnyValidation.sv_load_params(REPO_ROOT, controls)
    params = _apply_param_overrides(params_loaded.params)
    params_path = params_loaded.path
    param_overrides = _coerce_param_overrides(INTERACTIVE_PARAM_OVERRIDES[])

    fields = Float64.(_workflow_array(diag, "fields_T", controls["common"]["fields_T"]))
    Eis = Float64.(_workflow_array(diag, "Ei_meV", [get(controls["kpm"], "Ei_meV", 4.65)]))
    qtags = String.(_workflow_array(diag, "qtags", controls["kpm"]["qtags"]))
    include_flat = _as_bool(_workflow_get(diag, "include_flat_component", false))
    include_disp = _as_bool(_workflow_get(diag, "include_dispersive_component", true))
    include_disp || error("This first workflow expects include_dispersive_component=true")

    table_dir, fig_dir = _output_dirs(REPO_ROOT, diag)
    timers = TimerTable()

    println("Sunny KPM fixed-model summary")
    println("-----------------------------")
    println("repo_root      = ", REPO_ROOT)
    println("controls       = ", CONTROL_PATH)
    println("params         = ", params_path)
    println("Sunny version  = ", SunnyValidation.sv_try_pkgversion(Sunny))
    println("Ei_meV         = ", Eis)
    println("fields_T       = ", fields)
    println("qtags          = ", qtags)
    println("include_flat   = ", include_flat)
    println("scale_mode     = ", _workflow_get(diag, "neutron_scale_mode", "panel_least_squares"))
    println("KPM energy     = ", controls["kpm"]["energy_min_meV"], " to ",
            controls["kpm"]["energy_max_meV"], " meV, n = ", controls["kpm"]["n_energy"])
    println("KPM kernel     = ", controls["kpm"]["kernel_fwhm_meV"], " meV FWHM")
    println("E resample     = ", get(controls["kpm"], "energy_resample_mode", "repo_resolved"))
    if param_overrides !== nothing
        println("param_overrides = ", param_overrides)
    end

    backend = _maybe_make_backend(diag)
    println("KPM backend    = ", backend === nothing ? "CPU" : string(typeof(backend)))
    println()

    # KPM contexts are independent of Ei and qtag. Reuse one context per field/component.
    disp_contexts = Dict{Float64,Any}()
    flat_contexts = Dict{Float64,Any}()

    for B in fields
        disp_contexts[B] = _prepare_kpm_context(params, controls, diag, timers;
                                                component=:dispersive, field_T=B,
                                                backend=backend)
        if include_flat
            flat_contexts[B] = _prepare_kpm_context(params, controls, diag, timers;
                                                    component=:flat, field_T=B,
                                                    backend=backend)
        end
    end

    neutron_items = Any[]

    for Ei in Eis, B in fields, qtag in qtags
        cut, cut_controls = _try_load_cut(REPO_ROOT, controls; Ei_meV=Ei, field_T=B, qtag=qtag)
        cut === nothing && continue

        @info "Computing fixed-model KPM cut" Ei_meV=Ei field_T=B qtag=qtag
        disp = _component_cut_spectrum_from_context(disp_contexts[B], params, cut_controls, cut, timers)

        model_unscaled = copy(disp.intensity_model_expgrid)
        model_raw_unscaled = copy(disp.intensity_model_grid)
        component_mode = "dispersive"

        if include_flat
            flat = _component_cut_spectrum_from_context(flat_contexts[B], params, cut_controls, cut, timers)
            flat_weight = SunnyValidation.sv_flat_neutron_weight(params, cut_controls)
            model_unscaled .+= flat_weight .* flat.intensity_model_expgrid
            model_raw_unscaled .+= flat_weight .* flat.intensity_model_grid
            component_mode = "dispersive_plus_flat"
        end

        scale = _panel_scale(model_unscaled, cut.intensity, params, cut_controls, diag)
        model_scaled = scale .* model_unscaled
        model_raw_scaled = scale .* model_raw_unscaled

        push!(neutron_items, (;
            cut=cut,
            energy_meV=cut.energy_meV,
            intensity_exp=cut.intensity,
            model_unscaled,
            model_scaled,
            energy_model_meV=disp.energy_model_meV,
            model_raw_unscaled,
            model_raw_scaled,
            panel_scale=scale,
            q_samples=disp.q_samples,
            q_measured_samples=disp.q_measured_samples,
            q_resolution_samples=disp.q_resolution_samples,
            gpu_used=disp.gpu_used,
            component_mode,
            energy_resample_mode=String(get(cut_controls["kpm"], "energy_resample_mode", "repo_resolved")),
            form_factor_weight=disp.form_factor_weight,
            form_factor_amplitude=disp.form_factor_amplitude,
            qmag_Ainv=disp.qmag_Ainv,
        ))
    end

    # Mean-field bridge magnetization. This writes its own diagnostic CSV/PNG too,
    # but we also reuse the returned arrays in the combined summary figure.
    mag_result = SunnyValidation.sv_run_meanfield_magnetization(REPO_ROOT; controls=controls)

    # Keep the bulk panel consistent with the neutron component choice.
    # The helper above computes the canonical mean-field bridge including the
    # flat component; for this diagnostic workflow we can zero the flat
    # contribution in the combined summary figure when the neutron calculation
    # is dispersive-only.
    if !include_flat
        mag_result = merge(mag_result, (;
            M_total = mag_result.M_disp .+ mag_result.M_vv,
            M_flat = zeros(length(mag_result.B_T)),
            r2 = 0.0,
        ))
    end

    mag_path = SunnyValidation.sv_repo_path(REPO_ROOT, controls["paths"]["magnetization_csv"])
    mag_data = SunnyValidation.sv_read_magnetization_csv(mag_path)

    neutron_csv = joinpath(table_dir, "sunny_kpm_gpu_fixed_model_neutron_1d.csv")
    raw_neutron_csv = joinpath(table_dir, "sunny_kpm_gpu_fixed_model_neutron_1d_raw_energy_grid.csv")
    if !isempty(neutron_items)
        _write_neutron_long_csv(neutron_csv, neutron_items)
        _write_neutron_raw_energy_csv(raw_neutron_csv, neutron_items)
    end

    profile_csv = joinpath(table_dir, "sunny_kpm_gpu_fixed_model_profile.csv")
    _write_rows_csv(profile_csv, timers.rows)

    fig_path = joinpath(fig_dir, "sunny_kpm_gpu_fixed_model_summary.png")
    _plot_summary(fig_path, neutron_items, mag_result, mag_data, Eis, fields, qtags, diag)

    println()
    println("Saved:")
    println("  ", neutron_csv)
    println("  ", raw_neutron_csv)
    println("  ", profile_csv)
    println("  ", fig_path)
    println()
    println("Done.")
    return (; neutron_items, mag_result, neutron_csv, raw_neutron_csv, profile_csv, fig_path)
end

"""
    run_sunny_kpm_gpu_fixed_model_summary(; enable_gpu=nothing,
                                           backend=nothing,
                                           batched=nothing,
                                           precision=nothing,
                                           Ei_meV=nothing,
                                           neutron_intensity_ylim=nothing,
                                           include_flat_component=nothing,
                                           param_overrides=nothing,
                                           kpm_energy_min_meV=nothing,
                                           kpm_energy_max_meV=nothing,
                                           kpm_n_energy=nothing,
                                           kpm_kernel_fwhm_meV=nothing,
                                           energy_resample_mode=nothing)

Interactive-session wrapper for `main()`.

This is the recommended entry point when running from a warm Julia session, e.g.

```julia
include("scripts/sunny_kpm_gpu_fixed_model_summary.jl")
run_sunny_kpm_gpu_fixed_model_summary(; enable_gpu=true,
                                      backend="CUDA",
                                      batched=true,
                                      precision="Float64",
                                      Ei_meV=[4.65],
                                      neutron_intensity_ylim=(0.0, 0.003),
                                      include_flat_component=false,
                                      kpm_energy_max_meV=4.2,
                                      kpm_n_energy=241,
                                      energy_resample_mode="direct_interpolation",
                                      param_overrides=(; J1_meV=0.25,
                                                        J2_meV=0.01,
                                                        sigma_J=0.30,
                                                        gzz=3.8,
                                                        sigma_gzz=0.4))
```

Keyword arguments update the same environment variables used by the command-line
workflow before calling `main()`. Omitted keywords leave the current environment
unchanged.
"""
function run_sunny_kpm_gpu_fixed_model_summary(; enable_gpu=nothing,
                                                backend=nothing,
                                                batched=nothing,
                                                precision=nothing,
                                                Ei_meV=nothing,
                                                neutron_intensity_ylim=nothing,
                                                include_flat_component=nothing,
                                                param_overrides=nothing,
                                                kpm_energy_min_meV=nothing,
                                                kpm_energy_max_meV=nothing,
                                                kpm_n_energy=nothing,
                                                kpm_kernel_fwhm_meV=nothing,
                                                energy_resample_mode=nothing)
    if enable_gpu !== nothing
        ENV["SUNNY_KPM_ENABLE_GPU"] = enable_gpu ? "1" : "0"
    end
    if backend !== nothing
        ENV["SUNNY_KPM_GPU_BACKEND"] = String(backend)
    end
    if batched !== nothing
        ENV["SUNNY_KPM_GPU_BATCHED"] = batched ? "1" : "0"
    end
    if precision !== nothing
        ENV["SUNNY_KPM_GPU_PRECISION"] = String(precision)
    end

    overrides = Dict{String,Any}()
    if Ei_meV !== nothing
        overrides["Ei_meV"] = collect(Float64.(Ei_meV))
    end
    if neutron_intensity_ylim !== nothing
        overrides["neutron_intensity_ylim"] = collect(Float64.(neutron_intensity_ylim))
    end
    if include_flat_component !== nothing
        overrides["include_flat_component"] = Bool(include_flat_component)
    end
    if kpm_energy_min_meV !== nothing
        overrides["kpm_energy_min_meV"] = Float64(kpm_energy_min_meV)
    end
    if kpm_energy_max_meV !== nothing
        overrides["kpm_energy_max_meV"] = Float64(kpm_energy_max_meV)
    end
    if kpm_n_energy !== nothing
        overrides["kpm_n_energy"] = Int(kpm_n_energy)
    end
    if kpm_kernel_fwhm_meV !== nothing
        overrides["kpm_kernel_fwhm_meV"] = Float64(kpm_kernel_fwhm_meV)
    end
    if energy_resample_mode !== nothing
        overrides["energy_resample_mode"] = String(energy_resample_mode)
    end

    old_overrides = INTERACTIVE_WORKFLOW_OVERRIDES[]
    old_param_overrides = INTERACTIVE_PARAM_OVERRIDES[]
    INTERACTIVE_WORKFLOW_OVERRIDES[] = overrides
    INTERACTIVE_PARAM_OVERRIDES[] = param_overrides
    try
        return main()
    finally
        INTERACTIVE_WORKFLOW_OVERRIDES[] = old_overrides
        INTERACTIVE_PARAM_OVERRIDES[] = old_param_overrides
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
