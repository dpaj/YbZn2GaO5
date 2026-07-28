#!/usr/bin/env julia

# Sunny KPM diagnostic: dispersive component only, deterministic 1D Q-grid
# histogramming, and independent multipliers applied to the canonical best-fit
# exchange disorder (sigma_J) and dispersive g-factor disorder (sigma_gzz).
#
# Default diagnostic in configs/sunny_kpm_1d_disp_grid_2sigmaJ_controls.toml:
#   measured grid:  5 × 5 × 1 over the analytical 1D cut volume
#   resolution grid: 5 × 5 × 1 Gaussian momentum-resolution offsets
#   sigma_J used:   2 × best-fit sigma_J
#   sigma_gzz used: 1 × best-fit sigma_gzz
#
# This is a targeted script for testing whether the Sunny random-bond dispersive
# line shape can be reconciled with the analytical parameter scale, without the
# nondispersive/flat component.

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

# -----------------------------------------------------------------------------
# Small local helpers
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

function _load_diagnostic_controls(repo_root::AbstractString)
    diag_path = joinpath(repo_root, "configs", "sunny_kpm_1d_disp_grid_2sigmaJ_controls.toml")
    isfile(diag_path) || error("Could not find diagnostic controls: $diag_path")
    diag = load_toml_config(diag_path)

    base_rel = get(get(diag, "paths", Dict{String,Any}()), "base_controls_toml", "configs/sunny_validation_controls.toml")
    base_path = _repo_path(repo_root, String(base_rel))
    controls = load_toml_config(base_path)

    # Apply only the nested [control_overrides] table to the Sunny controls.
    if haskey(diag, "control_overrides")
        _deepmerge!(controls, diag["control_overrides"])
    end

    # Redirect output directories for this diagnostic while keeping the normal
    # data/input paths from the base Sunny controls.
    if haskey(diag, "paths")
        paths = diag["paths"]
        if haskey(paths, "figure_subdir")
            controls["paths"]["figure_subdir"] = paths["figure_subdir"]
        end
        if haskey(paths, "table_subdir")
            controls["paths"]["table_subdir"] = paths["table_subdir"]
        end
    end

    return (; diag, controls, diag_path, base_path)
end

function _replace_namedtuple(nt::NamedTuple; kwargs...)
    return merge(nt, (; kwargs...))
end


_ntget(nt, key::Symbol, default) = hasproperty(nt, key) ? getproperty(nt, key) : default

function _scale_mask(E, y, x, err; energy_window=(0.5, 3.0), positive_experiment=false, positive_model=false)
    mask = isfinite.(E) .& isfinite.(y) .& isfinite.(x) .& isfinite.(err)
    if energy_window !== nothing
        lo, hi = energy_window
        mask .&= (E .>= lo) .& (E .<= hi)
    end
    positive_experiment && (mask .&= y .> 0)
    positive_model && (mask .&= x .> 0)
    return mask
end

function _best_positive_scale(y, x, err, E;
        energy_window=(0.5, 3.0), use_uncertainties=true,
        nonnegative=true, positive_experiment=false, positive_model=false,
        fallback=1.0)
    mask = _scale_mask(E, y, x, err;
        energy_window, positive_experiment, positive_model)
    if !any(mask)
        @warn "No points available for scale fit; using fallback" fallback
        return Float64(fallback)
    end

    e = Float64.(err)
    if use_uncertainties
        finite_positive = e[isfinite.(e) .& (e .> 0)]
        floor = isempty(finite_positive) ? 1.0 : minimum(finite_positive)
        e = map(v -> (isfinite(v) && v > 0) ? v : floor, e)
        w = 1.0 ./ (e .^ 2)
    else
        w = ones(Float64, length(e))
    end

    xx = Float64.(x[mask])
    yy = Float64.(y[mask])
    ww = Float64.(w[mask])
    denom = sum(ww .* xx .* xx)
    if !(isfinite(denom) && denom > 0)
        @warn "Degenerate model for scale fit; using fallback" fallback denom
        return Float64(fallback)
    end
    s = sum(ww .* xx .* yy) / denom
    nonnegative && (s = max(0.0, s))
    return isfinite(s) ? Float64(s) : Float64(fallback)
end

function _scale_from_mode(scale_mode::Symbol, cut_rows, params, controls, diag)
    run = get(diag, "run", Dict{String,Any}())
    scale_fit_window = Tuple(Float64.(get(run, "scale_fit_window_meV", [0.5, 3.0])))
    use_uncertainties = Bool(get(run, "scale_use_uncertainties", true))
    positive_experiment = Bool(get(run, "scale_positive_experiment_only", false))
    positive_model = Bool(get(run, "scale_positive_model_only", false))
    manual_scale = Float64(get(run, "manual_neutron_scale", 1.0))
    best_fit_scale = SunnyValidation.sv_neutron_scale(params, controls)

    if scale_mode == :best_fit
        return Dict((r.field_T, r.qtag) => best_fit_scale for r in cut_rows), best_fit_scale
    elseif scale_mode == :manual
        return Dict((r.field_T, r.qtag) => manual_scale for r in cut_rows), manual_scale
    elseif scale_mode == :global_least_squares
        E = vcat([r.energy_meV for r in cut_rows]...)
        y = vcat([r.I_exp for r in cut_rows]...)
        x = vcat([r.I_model_unscaled for r in cut_rows]...)
        err = vcat([r.Ierr_exp for r in cut_rows]...)
        s = _best_positive_scale(y, x, err, E;
            energy_window=scale_fit_window,
            use_uncertainties,
            nonnegative=true,
            positive_experiment,
            positive_model,
            fallback=best_fit_scale)
        return Dict((r.field_T, r.qtag) => s for r in cut_rows), s
    elseif scale_mode == :panel_least_squares || scale_mode == :per_cut_least_squares
        d = Dict{Tuple{Float64,String},Float64}()
        for r in cut_rows
            d[(r.field_T, r.qtag)] = _best_positive_scale(r.I_exp, r.I_model_unscaled, r.Ierr_exp, r.energy_meV;
                energy_window=scale_fit_window,
                use_uncertainties,
                nonnegative=true,
                positive_experiment,
                positive_model,
                fallback=best_fit_scale)
        end
        return d, NaN
    else
        error("Unknown [run].scale_mode=$scale_mode. Use best_fit, manual, global_least_squares, or panel_least_squares.")
    end
end

function _qtag_index_map(qtags)
    return Dict(String(q) => i for (i, q) in enumerate(qtags))
end

function _write_results_csv(path, rows, scale_by_cut, best_fit_scale, params_orig, params_used, diag, controls)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join([
            "field_T", "qtag", "energy_meV",
            "I_exp", "Ierr_exp", "I_disp_scaled", "I_disp_unscaled", "I_disp_bestfit_scale", "residual",
            "neutron_scale", "best_fit_neutron_scale", "scale_mode",
            "sigma_J_original", "sigma_J_multiplier", "sigma_J_used",
            "sigma_gzz_original", "sigma_gzz_multiplier", "sigma_gzz_used",
            "q_samples", "q_measured_samples", "q_resolution_samples",
            "measured_n_h", "measured_n_k", "measured_n_l",
            "resolution_n_h", "resolution_n_k", "resolution_n_l",
            "sigma_H_rlu", "sigma_K_rlu", "sigma_L_rlu", "grid_nsigma",
            "include_exchange_disorder", "include_gzz_disorder",
            "sunny_dims_x", "sunny_dims_y", "sunny_dims_z",
            "sunny_repeat_x", "sunny_repeat_y", "sunny_repeat_z",
            "sunny_system_size_x", "sunny_system_size_y", "sunny_system_size_z",
        ], ","))

        run = diag["run"]
        scale_mode = String(Symbol(get(run, "scale_mode", "panel_least_squares")))
        sigmaJ_mult = Float64(get(run, "sigma_J_multiplier", 2.0))
        sigmag_mult = Float64(get(run, "sigma_gzz_multiplier", 1.0))
        sizectl = SunnyValidation.sv_system_size_controls(controls, "kpm")
        kc = controls["kpm"]

        for r in rows
            s = scale_by_cut[(r.field_T, r.qtag)]
            for i in eachindex(r.energy_meV)
                model_scaled = s * r.I_model_unscaled[i]
                model_bestfit = best_fit_scale * r.I_model_unscaled[i]
                residual = r.I_exp[i] - model_scaled
                vals = Any[
                    r.field_T, r.qtag, r.energy_meV[i],
                    r.I_exp[i], r.Ierr_exp[i], model_scaled, r.I_model_unscaled[i], model_bestfit, residual,
                    s, best_fit_scale, scale_mode,
                    params_orig.sigma_J, sigmaJ_mult, params_used.sigma_J,
                    params_orig.sigma_gzz, sigmag_mult, params_used.sigma_gzz,
                    r.q_samples, r.q_measured_samples, r.q_resolution_samples,
                    r.measured_n_h, r.measured_n_k, r.measured_n_l,
                    r.resolution_n_h, r.resolution_n_k, r.resolution_n_l,
                    r.sigma_H, r.sigma_K, r.sigma_L, r.grid_nsigma,
                    get(kc, "include_exchange_disorder", true), get(kc, "include_gzz_disorder", true),
                    sizectl.dims[1], sizectl.dims[2], sizectl.dims[3],
                    sizectl.repeat_factor[1], sizectl.repeat_factor[2], sizectl.repeat_factor[3],
                    sizectl.system_size[1], sizectl.system_size[2], sizectl.system_size[3],
                ]
                println(io, join(vals, ","))
            end
        end
    end
end

function _write_manifest_csv(path, rows, scale_by_cut, best_fit_scale, params_orig, params_used, diag, controls)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "field_T,qtag,n_energy,E_min,E_max,q_samples,q_measured_samples,q_resolution_samples,neutron_scale,best_fit_neutron_scale,sigma_J_original,sigma_J_multiplier,sigma_J_used,sigma_gzz_original,sigma_gzz_multiplier,sigma_gzz_used,include_exchange_disorder,include_gzz_disorder")
        run = diag["run"]
        sigmaJ_mult = Float64(get(run, "sigma_J_multiplier", 2.0))
        sigmag_mult = Float64(get(run, "sigma_gzz_multiplier", 1.0))
        for r in rows
            s = scale_by_cut[(r.field_T, r.qtag)]
            println(io, join(Any[
                r.field_T, r.qtag, length(r.energy_meV), minimum(r.energy_meV), maximum(r.energy_meV),
                r.q_samples, r.q_measured_samples, r.q_resolution_samples,
                s, best_fit_scale, params_orig.sigma_J, sigmaJ_mult, params_used.sigma_J,
                params_orig.sigma_gzz, sigmag_mult, params_used.sigma_gzz,
                get(controls["kpm"], "include_exchange_disorder", true),
                get(controls["kpm"], "include_gzz_disorder", true),
            ], ","))
        end
    end
end

function _make_plot(path, rows, scale_by_cut, best_fit_scale, diag, controls)
    mkpath(dirname(path))
    fields = Float64.(get(diag["run"], "fields_T", controls["common"]["fields_T"]))
    qtags = String.(get(diag["run"], "qtags", controls["kpm"]["qtags"]))
    qidx = _qtag_index_map(qtags)
    fidx = Dict(Float64(B) => i for (i, B) in enumerate(fields))

    fig = Figure(size=(430 * length(fields), 285 * length(qtags)))
    ylims = get(diag["run"], "plot_ylim", get(controls["kpm"], "plot_ylim", nothing))
    run = diag["run"]
    sigmaJ_mult = Float64(get(run, "sigma_J_multiplier", 2.0))
    sigmag_mult = Float64(get(run, "sigma_gzz_multiplier", 1.0))
    default_title = @sprintf("Sunny KPM dispersive only: 5×5 measured × 5×5 resolution, %.3g×σJ, %.3g×σgzz", sigmaJ_mult, sigmag_mult)
    title_suffix = String(get(run, "title_suffix", default_title))
    plot_best_fit_scale = Bool(get(run, "plot_best_fit_scale", false))

    for r in rows
        iq = qidx[r.qtag]
        iB = fidx[r.field_T]
        ax = Axis(fig[iq, iB], xlabel="Energy transfer (meV)", ylabel="Intensity", title=@sprintf("%s, %.0f T", r.qtag, r.field_T))
        scatter!(ax, r.energy_meV, r.I_exp; markersize=5, label="experiment")
        s = scale_by_cut[(r.field_T, r.qtag)]
        lines!(ax, r.energy_meV, s .* r.I_model_unscaled; label="Sunny disp., fitted scale")
        if plot_best_fit_scale
            lines!(ax, r.energy_meV, best_fit_scale .* r.I_model_unscaled; linestyle=:dash, label="Sunny disp., best-fit scale")
        end
        if ylims !== nothing && length(ylims) == 2
            ylims!(ax, Float64(ylims[1]), Float64(ylims[2]))
        end
        iq == 1 && iB == length(fields) && axislegend(ax, position=:rt)
    end
    Label(fig[0, :], title_suffix, fontsize=18)
    save(path, fig)
    return path
end


# -----------------------------------------------------------------------------
# Reusable Sunny/KPM context helpers for this diagnostic
# -----------------------------------------------------------------------------

function _timed!(profile_rows::Vector{NamedTuple}, field_T, qtag::AbstractString, stage::AbstractString, f;
        q_samples::Int=-1, n_energy::Int=-1, note::AbstractString="")
    local val
    t = @elapsed begin
        val = f()
    end
    push!(profile_rows, (;
        field_T = Float64(field_T),
        qtag = String(qtag),
        stage = String(stage),
        seconds = Float64(t),
        q_samples = Int(q_samples),
        n_energy = Int(n_energy),
        note = String(note),
    ))
    return val
end

function _write_profile_csv(path::AbstractString, profile_rows::Vector{NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "field_T,qtag,stage,seconds,q_samples,n_energy,note")
        for r in profile_rows
            println(io, join(Any[
                r.field_T, r.qtag, r.stage, @sprintf("%.9g", r.seconds),
                r.q_samples, r.n_energy, replace(r.note, "," => ";"),
            ], ","))
        end
    end
    return path
end

function _try_energy_per_site(sys)
    try
        return Float64(energy_per_site(sys))
    catch
        return NaN
    end
end

function _initialize_field_polarized!(sys, controls::Dict; field_T::Real=1.0)
    u = SunnyValidation.sv_field_direction(controls)
    sgn = Float64(field_T) < 0 ? -1.0 : 1.0
    dip = collect(sgn .* u)
    for site in eachsite(sys)
        set_dipole!(sys, dip, site)
    end
    return sys
end

function _build_reusable_dispersive_kpm_context(params, controls::Dict, field_T::Real, profile_rows::Vector{NamedTuple}, diag::Dict)
    run = get(diag, "run", Dict{String,Any}())
    kc = controls["kpm"]
    sizectl = SunnyValidation.sv_system_size_controls(controls, "kpm")
    include_exchange = Bool(get(kc, "include_exchange_disorder", true))
    include_gzz = Bool(get(kc, "include_gzz_disorder", true))
    maxiters = Int(kc["maxiters"])
    relax_ground_state = Bool(get(run, "relax_ground_state", true))
    initial_spin_state = Symbol(get(run, "initial_spin_state", "field_polarized"))

    state = _timed!(profile_rows, field_T, "ALL", "build_system", () -> begin
        base = SunnyValidation.sv_build_effective_sunny_system(params, controls;
            component=:dispersive, dims=sizectl.dims, field_T=field_T)
        sys = base.sys
        if sizectl.repeat_factor != (1, 1, 1)
            sys = to_inhomogeneous(repeat_periodically(sys, sizectl.repeat_factor))
        else
            sys = to_inhomogeneous(sys)
        end
        (; sys, crystal=base.crystal, units=base.units)
    end; note=@sprintf("dims=%s repeat=%s system_size=%s", string(sizectl.dims), string(sizectl.repeat_factor), string(sizectl.system_size)))

    sys = state.sys

    _timed!(profile_rows, field_T, "ALL", "apply_disorder", () -> begin
        SunnyValidation.sv_apply_disorder!(sys, params, controls;
            component=:dispersive, include_exchange=include_exchange, include_gzz=include_gzz)
        nothing
    end; note=@sprintf("include_exchange=%s include_gzz=%s sigma_J=%.6g sigma_gzz=%.6g", include_exchange, include_gzz, params.sigma_J, params.sigma_gzz))

    E_initial = _timed!(profile_rows, field_T, "ALL", "initialize_spins", () -> begin
        if initial_spin_state == :field_polarized
            _initialize_field_polarized!(sys, controls; field_T=field_T)
        elseif initial_spin_state == :random
            Random.seed!(Int(controls["common"]["seed"]) + 101)
            randomize_spins!(sys)
        elseif initial_spin_state in (:none, :as_built)
            # Keep Sunny's construction default.
        else
            error("Unknown [run].initial_spin_state=$(initial_spin_state). Use field_polarized, random, or none.")
        end
        _try_energy_per_site(sys)
    end; note=String(initial_spin_state))

    E_final = E_initial
    if relax_ground_state
        E_final = _timed!(profile_rows, field_T, "ALL", "minimize_energy", () -> begin
            minimize_energy!(sys; maxiters=maxiters)
            _try_energy_per_site(sys)
        end; note=@sprintf("maxiters=%d E_initial=%.10g", maxiters, E_initial))
    else
        push!(profile_rows, (;
            field_T = Float64(field_T), qtag = "ALL", stage = "minimize_energy",
            seconds = 0.0, q_samples = -1, n_energy = -1,
            note = @sprintf("skipped relax_ground_state=false E_initial=%.10g", E_initial),
        ))
    end

    energies = collect(range(Float64(kc["energy_min_meV"]), Float64(kc["energy_max_meV"]); length=Int(kc["n_energy"])))
    kernel = gaussian(fwhm=Float64(kc["kernel_fwhm_meV"]))

    swt = _timed!(profile_rows, field_T, "ALL", "construct_kpm", () -> begin
        measure = SunnyValidation.sv_sunny_measure(sys, controls)
        SpinWaveTheoryKPM(sys; measure=measure, tol=Float64(kc["tol"]))
    end; n_energy=length(energies), note=@sprintf("E_final=%.10g kernel_fwhm=%.6g tol=%.6g", E_final, Float64(kc["kernel_fwhm_meV"]), Float64(kc["tol"])))

    return (;
        sys = sys,
        swt = swt,
        energies = energies,
        kernel = kernel,
        field_T = Float64(field_T),
        component = :dispersive,
        initial_spin_state = initial_spin_state,
        relax_ground_state = relax_ground_state,
        E_initial = E_initial,
        E_final = E_final,
        sizectl = sizectl,
    )
end

function _kpm_component_spectrum_from_context(ctx, controls::Dict, cut, profile_rows::Vector{NamedTuple}; qtag::AbstractString)
    sampler = SunnyValidation.sv_kpm_1d_q_sampler(cut, controls)
    qs = sampler.qs
    energies = ctx.energies
    res = _timed!(profile_rows, ctx.field_T, qtag, "intensities", () -> begin
        intensities(ctx.swt, qs; energies=energies, kernel=ctx.kernel)
    end; q_samples=length(qs), n_energy=length(energies))

    raw = SunnyValidation.sv_try_extract_sunny_intensity(res)
    I0 = SunnyValidation.sv_orient_sunny_intensity_matrix(raw, length(energies), length(qs))
    form = _timed!(profile_rows, ctx.field_T, qtag, "form_factor_and_q_average", () -> begin
        Ipost, form_factor_weight, form_factor_amplitude, qmag_Ainv =
            SunnyValidation.sv_apply_form_factor_to_intensity(I0, qs, controls)
        Iavg = SunnyValidation.sv_kpm_1d_average_qsampled_intensity(Ipost, sampler)
        (; Ipost, I0, Iavg, form_factor_weight, form_factor_amplitude, qmag_Ainv)
    end; q_samples=length(qs), n_energy=length(energies))

    return (;
        energy_meV = energies,
        intensity = form.Iavg,
        intensity_qsampled = form.Ipost,
        intensity_no_form_factor = form.I0,
        result = res,
        q = Float64.(sampler.q_center),
        qs = qs,
        q_average = sampler,
        q_average_enabled = sampler.q_average_enabled,
        q_samples = sampler.n_samples,
        q_measured_samples = sampler.n_measured,
        q_resolution_samples = sampler.n_resolution,
        form_factor_weight = sum(sampler.weights .* form.form_factor_weight),
        form_factor_amplitude = sum(sampler.weights .* form.form_factor_amplitude),
        qmag_Ainv = sum(sampler.weights .* form.qmag_Ainv),
    )
end

# -----------------------------------------------------------------------------
# Main run
# -----------------------------------------------------------------------------

function main()
    (; diag, controls, diag_path, base_path) = _load_diagnostic_controls(REPO_ROOT)
    params_loaded = SunnyValidation.sv_load_params(REPO_ROOT, controls)
    params = params_loaded.params
    params_path = params_loaded.path

    run = get(diag, "run", Dict{String,Any}())
    sigmaJ_mult = Float64(get(run, "sigma_J_multiplier", 2.0))
    sigmag_mult = Float64(get(run, "sigma_gzz_multiplier", 1.0))
    params2 = _replace_namedtuple(params;
        sigma_J = sigmaJ_mult * params.sigma_J,
        sigma_gzz = sigmag_mult * params.sigma_gzz)

    fields = Float64.(get(run, "fields_T", controls["common"]["fields_T"]))
    qtags = String.(get(run, "qtags", controls["kpm"]["qtags"]))
    scale_mode = Symbol(get(run, "scale_mode", "panel_least_squares"))
    hist_mode = Symbol(get(controls["kpm"], "histogram_mode", "bin_average"))

    out_table_dir = _repo_path(REPO_ROOT, controls["paths"]["table_subdir"])
    out_fig_dir = _repo_path(REPO_ROOT, controls["paths"]["figure_subdir"])
    mkpath(out_table_dir)
    mkpath(out_fig_dir)

    println("Sunny dispersive-only KPM diagnostic")
    println("------------------------------------")
    println("repo_root          = ", REPO_ROOT)
    println("controls           = ", diag_path)
    println("base controls      = ", base_path)
    println("params             = ", params_path)
    println("fields_T           = ", fields)
    println("qtags              = ", qtags)
    println("sigma_J original   = ", params.sigma_J)
    println("sigma_J multiplier = ", sigmaJ_mult)
    println("sigma_J used       = ", params2.sigma_J)
    println("sigma_gzz original = ", params.sigma_gzz)
    println("sigma_gzz mult.    = ", sigmag_mult)
    println("sigma_gzz used     = ", params2.sigma_gzz)
    println("scale_mode         = ", scale_mode)
    println()

    cuts_all = SunnyValidation.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
    cuts = [c for c in cuts_all if any(abs(c.field_T - B) < 1e-6 for B in fields) && c.qtag in qtags]
    isempty(cuts) && error("No experimental cuts matched fields=$fields qtags=$qtags")
    sort!(cuts; by = c -> (findfirst(==(c.qtag), qtags), findfirst(B -> abs(c.field_T - B) < 1e-6, fields)))

    rows = NamedTuple[]
    profile_rows = NamedTuple[]
    reuse_ground_state = Bool(get(run, "reuse_ground_state_by_field", true))

    if reuse_ground_state
        @info "Using one field-polarized/minimized Sunny ground state and one KPM object per field" fields sigma_J=params2.sigma_J sigma_gzz=params2.sigma_gzz
        context_by_field = Dict{Float64,Any}()
        for B in fields
            @info "Preparing reusable dispersive Sunny KPM context" field_T=B sigma_J=params2.sigma_J sigma_gzz=params2.sigma_gzz
            context_by_field[Float64(B)] = _build_reusable_dispersive_kpm_context(params2, controls, B, profile_rows, diag)
        end

        for cut in cuts
            Bkey = fields[findfirst(B -> abs(Float64(B) - cut.field_T) < 1e-6, fields)]
            ctx = context_by_field[Float64(Bkey)]
            @info "Computing dispersive Sunny KPM grid cut from reusable context" field_T=cut.field_T qtag=cut.qtag sigma_J=params2.sigma_J sigma_gzz=params2.sigma_gzz q_samples="grid"
            disp = _kpm_component_spectrum_from_context(ctx, controls, cut, profile_rows; qtag=cut.qtag)
            I_on_exp = SunnyValidation.sv_model_to_experimental_energy_grid_resolved(
                disp.energy_meV, disp.intensity, cut.energy_meV, controls; section="kpm", mode=hist_mode)

            qa = _ntget(disp, :q_average, (;))
            push!(rows, (;
                field_T = Float64(cut.field_T),
                qtag = String(cut.qtag),
                energy_meV = Float64.(cut.energy_meV),
                I_exp = Float64.(cut.intensity),
                Ierr_exp = Float64.(cut.error),
                I_model_unscaled = Float64.(I_on_exp),
                q_samples = Int(_ntget(disp, :q_samples, length(_ntget(disp, :qs, [])))),
                q_measured_samples = Int(_ntget(disp, :q_measured_samples, -1)),
                q_resolution_samples = Int(_ntget(disp, :q_resolution_samples, -1)),
                measured_n_h = Int(_ntget(qa, :measured_n_h, -1)),
                measured_n_k = Int(_ntget(qa, :measured_n_k, -1)),
                measured_n_l = Int(_ntget(qa, :measured_n_l, -1)),
                resolution_n_h = Int(_ntget(qa, :n_h, -1)),
                resolution_n_k = Int(_ntget(qa, :n_k, -1)),
                resolution_n_l = Int(_ntget(qa, :n_l, -1)),
                sigma_H = Float64(_ntget(qa, :sigma_H, NaN)),
                sigma_K = Float64(_ntget(qa, :sigma_K, NaN)),
                sigma_L = Float64(_ntget(qa, :sigma_L, NaN)),
                grid_nsigma = Float64(_ntget(qa, :grid_nsigma, NaN)),
            ))
        end
    else
        @info "Using legacy per-cut Sunny build/minimize path" sigma_J=params2.sigma_J sigma_gzz=params2.sigma_gzz
        for cut in cuts
            @info "Computing dispersive Sunny KPM grid cut" field_T=cut.field_T qtag=cut.qtag sigma_J=params2.sigma_J sigma_gzz=params2.sigma_gzz
            disp = _timed!(profile_rows, cut.field_T, cut.qtag, "legacy_sv_kpm_component_spectrum", () -> begin
                SunnyValidation.sv_kpm_component_spectrum(params2, controls;
                    component=:dispersive, field_T=cut.field_T, qtag=cut.qtag, cut=cut)
            end)
            I_on_exp = SunnyValidation.sv_model_to_experimental_energy_grid_resolved(
                disp.energy_meV, disp.intensity, cut.energy_meV, controls; section="kpm", mode=hist_mode)

            qa = _ntget(disp, :q_average, (;))
            push!(rows, (;
                field_T = Float64(cut.field_T),
                qtag = String(cut.qtag),
                energy_meV = Float64.(cut.energy_meV),
                I_exp = Float64.(cut.intensity),
                Ierr_exp = Float64.(cut.error),
                I_model_unscaled = Float64.(I_on_exp),
                q_samples = Int(_ntget(disp, :q_samples, length(_ntget(disp, :qs, [])))),
                q_measured_samples = Int(_ntget(disp, :q_measured_samples, -1)),
                q_resolution_samples = Int(_ntget(disp, :q_resolution_samples, -1)),
                measured_n_h = Int(_ntget(qa, :measured_n_h, -1)),
                measured_n_k = Int(_ntget(qa, :measured_n_k, -1)),
                measured_n_l = Int(_ntget(qa, :measured_n_l, -1)),
                resolution_n_h = Int(_ntget(qa, :n_h, -1)),
                resolution_n_k = Int(_ntget(qa, :n_k, -1)),
                resolution_n_l = Int(_ntget(qa, :n_l, -1)),
                sigma_H = Float64(_ntget(qa, :sigma_H, NaN)),
                sigma_K = Float64(_ntget(qa, :sigma_K, NaN)),
                sigma_L = Float64(_ntget(qa, :sigma_L, NaN)),
                grid_nsigma = Float64(_ntget(qa, :grid_nsigma, NaN)),
            ))
        end
    end

    scale_by_cut, global_scale = _scale_from_mode(scale_mode, rows, params, controls, diag)
    best_fit_scale = SunnyValidation.sv_neutron_scale(params, controls)

    csv_path = joinpath(out_table_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_vs_exp.csv")
    manifest_path = joinpath(out_table_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_manifest.csv")
    fig_path = joinpath(out_fig_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_vs_exp.png")
    profile_path = joinpath(out_table_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_profile.csv")

    _write_results_csv(csv_path, rows, scale_by_cut, best_fit_scale, params, params2, diag, controls)
    _write_manifest_csv(manifest_path, rows, scale_by_cut, best_fit_scale, params, params2, diag, controls)
    _make_plot(fig_path, rows, scale_by_cut, best_fit_scale, diag, controls)
    _write_profile_csv(profile_path, profile_rows)

    println()
    println("Wrote:")
    println("  ", csv_path)
    println("  ", manifest_path)
    println("  ", fig_path)
    println("  ", profile_path)
    println()
    println("Note: fitted-scale curve is for line-shape comparison. Set [run].plot_best_fit_scale=true to also draw the best-fit analytical scale.")

    return (; csv_path, manifest_path, fig_path, profile_path, rows, profile_rows, scale_by_cut, global_scale)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
