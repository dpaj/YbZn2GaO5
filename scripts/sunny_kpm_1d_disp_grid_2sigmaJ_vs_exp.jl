#!/usr/bin/env julia

# Sunny KPM diagnostic: dispersive component only, deterministic 1D Q-grid
# histogramming, and an exchange-disorder multiplier applied to the canonical
# best-fit sigma_J.
#
# Default diagnostic in configs/sunny_kpm_1d_disp_grid_2sigmaJ_controls.toml:
#   measured grid:  5 × 5 × 1 over the analytical 1D cut volume
#   resolution grid: 5 × 5 × 1 Gaussian momentum-resolution offsets
#   sigma_J used:   2 × best-fit sigma_J
#
# This is a targeted script for testing whether the Sunny random-bond dispersive
# line shape can be reconciled with the analytical parameter scale, without the
# nondispersive/flat component.

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using CairoMakie

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
        sigma_mult = Float64(get(run, "sigma_J_multiplier", 2.0))
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
                    params_orig.sigma_J, sigma_mult, params_used.sigma_J,
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
        println(io, "field_T,qtag,n_energy,E_min,E_max,q_samples,q_measured_samples,q_resolution_samples,neutron_scale,best_fit_neutron_scale,sigma_J_original,sigma_J_used,include_exchange_disorder,include_gzz_disorder")
        for r in rows
            s = scale_by_cut[(r.field_T, r.qtag)]
            println(io, join(Any[
                r.field_T, r.qtag, length(r.energy_meV), minimum(r.energy_meV), maximum(r.energy_meV),
                r.q_samples, r.q_measured_samples, r.q_resolution_samples,
                s, best_fit_scale, params_orig.sigma_J, params_used.sigma_J,
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
    title_suffix = String(get(diag["run"], "title_suffix", "Sunny dispersive only, 5×5 measured × 5×5 resolution, 2×σJ"))

    for r in rows
        iq = qidx[r.qtag]
        iB = fidx[r.field_T]
        ax = Axis(fig[iq, iB], xlabel="Energy transfer (meV)", ylabel="Intensity", title=@sprintf("%s, %.0f T", r.qtag, r.field_T))
        scatter!(ax, r.energy_meV, r.I_exp; markersize=5, label="experiment")
        s = scale_by_cut[(r.field_T, r.qtag)]
        lines!(ax, r.energy_meV, s .* r.I_model_unscaled; label="Sunny disp., fitted scale")
        lines!(ax, r.energy_meV, best_fit_scale .* r.I_model_unscaled; linestyle=:dash, label="Sunny disp., best-fit scale")
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
# Main run
# -----------------------------------------------------------------------------

function main()
    (; diag, controls, diag_path, base_path) = _load_diagnostic_controls(REPO_ROOT)
    params_loaded = SunnyValidation.sv_load_params(REPO_ROOT, controls)
    params = params_loaded.params
    params_path = params_loaded.path

    run = get(diag, "run", Dict{String,Any}())
    sigma_mult = Float64(get(run, "sigma_J_multiplier", 2.0))
    params2 = _replace_namedtuple(params; sigma_J = sigma_mult * params.sigma_J)

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
    println("sigma_J multiplier = ", sigma_mult)
    println("sigma_J used       = ", params2.sigma_J)
    println("scale_mode         = ", scale_mode)
    println()

    cuts_all = SunnyValidation.sv_load_kpm_experimental_cuts(REPO_ROOT, controls)
    cuts = [c for c in cuts_all if any(abs(c.field_T - B) < 1e-6 for B in fields) && c.qtag in qtags]
    isempty(cuts) && error("No experimental cuts matched fields=$fields qtags=$qtags")
    sort!(cuts; by = c -> (findfirst(==(c.qtag), qtags), findfirst(B -> abs(c.field_T - B) < 1e-6, fields)))

    rows = NamedTuple[]
    for cut in cuts
        @info "Computing dispersive Sunny KPM grid cut" field_T=cut.field_T qtag=cut.qtag sigma_J=params2.sigma_J
        disp = SunnyValidation.sv_kpm_component_spectrum(params2, controls;
            component=:dispersive, field_T=cut.field_T, qtag=cut.qtag, cut=cut)
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

    scale_by_cut, global_scale = _scale_from_mode(scale_mode, rows, params, controls, diag)
    best_fit_scale = SunnyValidation.sv_neutron_scale(params, controls)

    csv_path = joinpath(out_table_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_vs_exp.csv")
    manifest_path = joinpath(out_table_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_manifest.csv")
    fig_path = joinpath(out_fig_dir, "sunny_kpm_1d_disp_grid_2sigmaJ_vs_exp.png")

    _write_results_csv(csv_path, rows, scale_by_cut, best_fit_scale, params, params2, diag, controls)
    _write_manifest_csv(manifest_path, rows, scale_by_cut, best_fit_scale, params, params2, diag, controls)
    _make_plot(fig_path, rows, scale_by_cut, best_fit_scale, diag, controls)

    println()
    println("Wrote:")
    println("  ", csv_path)
    println("  ", manifest_path)
    println("  ", fig_path)
    println()
    println("Note: fitted-scale curve is for line-shape comparison. The dashed curve uses neutron_global_scale from best_fit_parameters.toml.")

    return (; csv_path, manifest_path, fig_path, rows, scale_by_cut, global_scale)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
