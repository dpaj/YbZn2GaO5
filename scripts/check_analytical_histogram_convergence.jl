# scripts/check_analytical_histogram_convergence.jl
#
# Analytical 1D histogram/event-sampling convergence study for YbZn2GaO5.
#
# This script does not optimize. It loads the canonical best-fit parameters,
# reruns the analytical neutron two-kernel model for a sequence of Monte Carlo
# sample counts, and compares each count against a high-sample reference.
#
# The aim is to estimate how many Q/disorder/event samples are needed before the
# experimental histogramming part of the analytical model is stable. That number
# can then guide expensive Sunny KPM histogramming tests.
#
# Run from the repo root:
#
#   julia --project=. scripts/check_analytical_histogram_convergence.jl

using Printf
using Statistics
using LinearAlgebra

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "scripts", "legacy", "YZGO_cofit_9T14T_shared_fraction_legacy.jl"))

function _get_nested(config::Dict, keys::Vector{String}, default=nothing)
    cur = config
    for (i, key) in enumerate(keys)
        if !(cur isa Dict) || !haskey(cur, key)
            return default
        end
        cur = cur[key]
    end
    return cur
end

function _as_string_vector(x)
    x === nothing && return String[]
    return String[string(v) for v in x]
end

function _as_float_vector(x)
    x === nothing && return Float64[]
    return Float64[Float64(v) for v in x]
end

function _as_int_vector(x)
    x === nothing && return Int[]
    return Int[round(Int, v) for v in x]
end

function _fit_param_dict_from_specs(specs)
    p = Dict{Symbol,Float64}()
    for s in specs
        p[s.name] = Float64(s.initial)
    end
    return p
end

function _model_key_for_qtag(qtag::String)
    return NF.OVERLAY_QTAG_TO_MODELKEY[qtag]
end

function _component_arrays(nres, field_T::Real, qtag::String; scale::Real=nres.scale, r2::Real=nres.r2)
    B = Float64(field_T)
    model_components = nres.model[B]
    model_key = _model_key_for_qtag(qtag)
    disp = model_components.dispersive[model_key]
    flat = model_components.nondispersive[model_key]
    E = Float64.(disp.E_centers_meV)
    Idisp = Float64(scale) .* Float64.(disp.intensity)
    Iflat = Float64(scale) .* Float64(r2) .* Float64.(flat.intensity)
    Itotal = Idisp .+ Iflat
    return (; E=E, total=Itotal, dispersive=Idisp, flat=Iflat,
            dispersive_unscaled=Float64.(disp.intensity), flat_unscaled=Float64.(flat.intensity))
end

function _rmse(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || error("Vector length mismatch: $(length(a)) vs $(length(b))")
    isempty(a) && return NaN
    return sqrt(mean((Float64.(a) .- Float64.(b)).^2))
end

function _relative_rmse(a::AbstractVector, b::AbstractVector)
    denom = sqrt(mean(Float64.(b).^2))
    return _rmse(a, b) / max(denom, eps(Float64))
end

function _integral(E::AbstractVector, I::AbstractVector)
    length(E) == length(I) || error("E/I length mismatch")
    length(E) < 2 && return NaN
    s = 0.0
    for i in 1:length(E)-1
        s += 0.5 * (I[i] + I[i+1]) * (E[i+1] - E[i])
    end
    return s
end

function _peak_metrics(E::AbstractVector, I::AbstractVector)
    if isempty(E)
        return (; peak_E=NaN, peak_I=NaN, centroid_E=NaN, linewidth_sigma=NaN, integrated=NaN)
    end
    If = Float64.(I)
    Ef = Float64.(E)
    imax = argmax(If)
    # Positive-weight centroid/linewidth proxy.  This is a robust diagnostic,
    # not a fit to a peak shape.
    W = max.(If, 0.0)
    wsum = sum(W)
    if wsum > 0
        mu = sum(W .* Ef) / wsum
        sig = sqrt(max(sum(W .* (Ef .- mu).^2) / wsum, 0.0))
    else
        mu = NaN
        sig = NaN
    end
    return (; peak_E=Ef[imax], peak_I=If[imax], centroid_E=mu, linewidth_sigma=sig,
            integrated=_integral(Ef, If))
end

function _write_csv(path::AbstractString, header::Vector{String}, rows::Vector{<:Tuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
    return path
end

function _write_manifest(path::AbstractString; rows)
    header = [
        "n_samples", "field_T", "qtag", "component",
        "peak_E_meV", "peak_I", "centroid_E_meV", "linewidth_sigma_meV", "integrated_intensity",
    ]
    _write_csv(path, header, rows)
end

function _write_metrics(path::AbstractString; rows)
    header = [
        "n_samples", "reference_n_samples", "field_T", "qtag", "component",
        "rmse", "relative_rmse", "max_abs_diff", "integrated", "reference_integrated",
        "relative_integrated_error", "peak_E_meV", "reference_peak_E_meV", "peak_E_shift_meV",
        "peak_I", "reference_peak_I", "relative_peak_height_error",
        "linewidth_sigma_meV", "reference_linewidth_sigma_meV", "relative_linewidth_error",
    ]
    _write_csv(path, header, rows)
end

function _write_model_long_csv(path::AbstractString, results_by_n::Dict{Int,Any}, fields_T, qtags)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "n_samples,field_T,qtag,energy_meV,I_total_scaled,I_dispersive_scaled,I_flat_scaled,I_dispersive_unscaled,I_flat_unscaled,neutron_global_scale,r2_shared")
        for n in sort(collect(keys(results_by_n)))
            nres = results_by_n[n]
            for B in fields_T, qtag in qtags
                arr = _component_arrays(nres, B, qtag)
                for i in eachindex(arr.E)
                    println(io, join((
                        n,
                        Float64(B),
                        qtag,
                        arr.E[i],
                        arr.total[i],
                        arr.dispersive[i],
                        arr.flat[i],
                        arr.dispersive_unscaled[i],
                        arr.flat_unscaled[i],
                        nres.scale,
                        nres.r2,
                    ), ","))
                end
            end
        end
    end
    return path
end

function _make_convergence_plots(figdir::AbstractString, results_by_n::Dict{Int,Any}, metrics_rows, fields_T, qtags, overlay_sample_counts)
    try
        @eval import CairoMakie as CM
    catch err
        @warn "Could not load CairoMakie; skipping convergence plots" exception=(err, catch_backtrace())
        return String[]
    end

    mkpath(figdir)
    files = String[]
    ns = sort(collect(keys(results_by_n)))
    ref_n = maximum(ns)

    # Overlay total 1D spectra for selected sample counts.
    overlay_ns = [n for n in overlay_sample_counts if haskey(results_by_n, n)]
    if !(ref_n in overlay_ns)
        push!(overlay_ns, ref_n)
    end
    overlay_ns = unique(overlay_ns)

    fig = CM.Figure(size=(1500, 800))
    for (ir, B) in enumerate(fields_T)
        for (ic, qtag) in enumerate(qtags)
            ax = CM.Axis(fig[ir, ic];
                title=@sprintf("%s, %.0f T", qtag, Float64(B)),
                xlabel="Energy transfer (meV)",
                ylabel=ic == 1 ? "Analytical model intensity" : "",
            )
            for n in overlay_ns
                arr = _component_arrays(results_by_n[n], B, qtag)
                lw = n == ref_n ? 3.0 : 1.5
                CM.lines!(ax, arr.E, arr.total; label="N=$(n)", linewidth=lw)
            end
            if ir == 1 && ic == length(qtags)
                CM.axislegend(ax; position=:rt, framevisible=false)
            end
        end
    end
    path_overlay = joinpath(figdir, "analytical_histogram_convergence_1d_overlays.png")
    CM.save(path_overlay, fig)
    push!(files, path_overlay)

    # Metric plot: relative RMSE of total component against reference.
    # Rebuild arrays from metric rows tuple layout.
    fig2 = CM.Figure(size=(1100, 700))
    ax2 = CM.Axis(fig2[1, 1];
        title="Analytical 1D histogram convergence vs reference",
        xlabel="log10(samples per cut per component)",
        ylabel="relative RMSE vs reference",
    )
    for B in fields_T, qtag in qtags
        xs = Float64[]
        ys = Float64[]
        for row in metrics_rows
            n = row[1]
            comp = row[5]
            if Float64(row[3]) == Float64(B) && row[4] == qtag && comp == "total" && n != ref_n
                push!(xs, log10(Float64(n)))
                push!(ys, Float64(row[7]))
            end
        end
        if !isempty(xs)
            CM.lines!(ax2, xs, ys; label=@sprintf("%s %.0fT", qtag, Float64(B)), linewidth=2)
            CM.scatter!(ax2, xs, ys; markersize=9)
        end
    end
    CM.axislegend(ax2; position=:rt, framevisible=false)
    path_metrics = joinpath(figdir, "analytical_histogram_convergence_relative_rmse.png")
    CM.save(path_metrics, fig2)
    push!(files, path_metrics)

    return files
end

function main()
    controls_path = joinpath(REPO_ROOT, "configs", "analytical_histogram_convergence_controls.toml")
    cofit_controls_path = joinpath(REPO_ROOT, "configs", "cofit_controls.toml")
    best_fit_path = joinpath(REPO_ROOT, "configs", "best_fit_parameters.toml")

    controls = load_toml_config(controls_path)
    cofit_controls = load_cofit_controls(cofit_controls_path)
    best_fit_config = load_best_fit_parameters(best_fit_path)

    initial_guess_kwargs = cofit_initial_guess_kwargs(best_fit_config)
    specs = cofit_default_param_specs(; initial_guess_kwargs...)
    pbest = _fit_param_dict_from_specs(specs)

    sample_counts = _as_int_vector(_get_nested(controls, ["sampling", "sample_counts"], [1000, 5000, 20000, 100000, 500000]))
    reference_n = round(Int, _get_nested(controls, ["sampling", "reference_n_samples"], maximum(sample_counts)))
    if !(reference_n in sample_counts)
        push!(sample_counts, reference_n)
    end
    sample_counts = sort(unique(sample_counts))

    seed = round(Int, _get_nested(controls, ["sampling", "seed"], 2026))
    table_dir = joinpath(REPO_ROOT, string(_get_nested(controls, ["output", "table_subdir"], "results/feature_tables/analytical_histogram_convergence")))
    fig_dir = joinpath(REPO_ROOT, string(_get_nested(controls, ["output", "figure_subdir"], "results/figures/analytical_histogram_convergence")))
    mkpath(table_dir)
    mkpath(fig_dir)

    fields_T = _as_float_vector(cofit_controls["data"]["fields_T"])
    qtags = _as_string_vector(cofit_controls["data"]["qtags"])

    model_controls = controls["model"]
    S = Float64(get(model_controls, "S", 0.5))
    J_units = toml_symbol(get(model_controls, "J_units", "fractional"))
    correlate_J1_J2 = Bool(get(model_controls, "correlate_J1_J2", false))
    use_form_factor = Bool(get(model_controls, "use_form_factor", true))
    include_j2_formfactor = Bool(get(model_controls, "include_j2_formfactor", true))
    include_kfki = Bool(get(model_controls, "include_kfki", true))
    polarization = toml_symbol(get(model_controls, "polarization", "transverse_c"))

    base_dir = joinpath(REPO_ROOT, "data", "neutron", "CNCS_1d_scans")
    data_mode = toml_symbol(cofit_controls["data"]["data_mode"])

    println("Analytical histogram convergence study")
    println("--------------------------------------")
    println("controls:        ", controls_path)
    println("best parameters: ", best_fit_path)
    println("data mode:       ", data_mode)
    println("sample counts:   ", sample_counts)
    println("reference N:     ", reference_n)
    println("fields:          ", fields_T)
    println("qtags:           ", qtags)
    println("output tables:   ", table_dir)
    println("output figures:  ", fig_dir)
    println()

    print_initial_guess_kwargs(initial_guess_kwargs)

    data_scans, _, _ = load_neutron_fit_data_1d_filtered(;
        base_dir=base_dir,
        data_mode=data_mode,
        Ei_meV=cofit_controls["data"]["neutron_fit_Ei_meV"],
        temperature_K=cofit_controls["data"]["neutron_fit_temperature_K"],
    )

    results_by_n = Dict{Int,Any}()

    for n in sample_counts
        println()
        println(@sprintf("Running analytical neutron model with n_samples_per_cut = %d", n))
        nres = evaluate_neutron_fit(pbest, data_scans;
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
            n_samples_per_cut=n,
            seed=seed,
            use_errors=true,
            error_floor=0.0,
            lattice=NF.demo_defaults().lattice,
            resolution=NF.demo_defaults().resolution,
            cuts=NF.default_cuts_1d(),
            S=S,
            J_units=J_units,
            correlate_J1_J2=correlate_J1_J2,
            use_form_factor=use_form_factor,
            include_j2_formfactor=include_j2_formfactor,
            include_kfki=include_kfki,
            polarization=polarization,
            nfree=0,
        )
        results_by_n[n] = nres

        NF.write_scaled_model_csv(joinpath(table_dir, @sprintf("analytical_1d_model_scaled_N%d.csv", n)), nres.model;
            scale=nres.scale,
            r2=nres.r2,
            qtags=qtags,
        )
        NF.write_fit_points_csv(joinpath(table_dir, @sprintf("analytical_1d_fit_points_N%d.csv", n)), nres.bundle, nres.scale, nres.r2)
        println(@sprintf("  redchi2 = %.8g, scale = %.8g, r2 = %.8g", nres.redchi2, nres.scale, nres.r2))
    end

    haskey(results_by_n, reference_n) || error("Reference sample count $reference_n was not run")
    ref = results_by_n[reference_n]

    manifest_rows = Tuple[]
    metrics_rows = Tuple[]
    components = _as_string_vector(_get_nested(controls, ["comparison", "components"], ["total", "dispersive", "flat"]))

    for n in sample_counts
        nres = results_by_n[n]
        for B in fields_T, qtag in qtags
            arr = _component_arrays(nres, B, qtag)
            arr_ref = _component_arrays(ref, B, qtag)
            for comp in components
                y = getproperty(arr, Symbol(comp))
                yref = getproperty(arr_ref, Symbol(comp))
                m = _peak_metrics(arr.E, y)
                mref = _peak_metrics(arr_ref.E, yref)
                push!(manifest_rows, (
                    n, Float64(B), qtag, comp,
                    m.peak_E, m.peak_I, m.centroid_E, m.linewidth_sigma, m.integrated,
                ))
                rmse = _rmse(y, yref)
                relrmse = _relative_rmse(y, yref)
                maxabs = maximum(abs.(Float64.(y) .- Float64.(yref)))
                relint = (m.integrated - mref.integrated) / max(abs(mref.integrated), eps(Float64))
                relpeak = (m.peak_I - mref.peak_I) / max(abs(mref.peak_I), eps(Float64))
                rellw = (m.linewidth_sigma - mref.linewidth_sigma) / max(abs(mref.linewidth_sigma), eps(Float64))
                push!(metrics_rows, (
                    n, reference_n, Float64(B), qtag, comp,
                    rmse, relrmse, maxabs,
                    m.integrated, mref.integrated, relint,
                    m.peak_E, mref.peak_E, m.peak_E - mref.peak_E,
                    m.peak_I, mref.peak_I, relpeak,
                    m.linewidth_sigma, mref.linewidth_sigma, rellw,
                ))
            end
        end
    end

    manifest_path = joinpath(table_dir, "analytical_histogram_convergence_peak_manifest.csv")
    metrics_path = joinpath(table_dir, "analytical_histogram_convergence_metrics_vs_reference.csv")
    long_path = joinpath(table_dir, "analytical_histogram_convergence_models_long.csv")
    _write_manifest(manifest_path; rows=manifest_rows)
    _write_metrics(metrics_path; rows=metrics_rows)
    _write_model_long_csv(long_path, results_by_n, fields_T, qtags)

    println()
    println("Wrote convergence tables:")
    println("  ", manifest_path)
    println("  ", metrics_path)
    println("  ", long_path)

    if Bool(_get_nested(controls, ["plotting", "make_plots"], true))
        overlay_counts = _as_int_vector(_get_nested(controls, ["plotting", "overlay_sample_counts"], [first(sample_counts), reference_n]))
        plot_files = _make_convergence_plots(fig_dir, results_by_n, metrics_rows, fields_T, qtags, overlay_counts)
        if !isempty(plot_files)
            println("Wrote convergence figures:")
            for f in plot_files
                println("  ", f)
            end
        end
    end

    println()
    println("Convergence study complete.")
    return (; results_by_n=results_by_n, metrics_path=metrics_path, manifest_path=manifest_path, long_path=long_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
