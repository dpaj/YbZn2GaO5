# scripts/check_analytical_histogram_grid_convergence.jl
#
# Deterministic-grid convergence study for the analytical 1D histogrammer.
#
# This is a companion to check_analytical_histogram_convergence.jl.  Instead of
# sampling measured Q and momentum resolution by Monte Carlo, it replaces those
# parts with deterministic H/K grids.  By default L is held at the center of the
# cut/resolution grid, because the field-polarized dispersion has no L
# dependence in this effective triangular-lattice model.
#
# The purpose is to test whether a modest deterministic Q grid can reproduce the
# high-statistics analytical Monte Carlo histogram.  If so, the same grid sizes
# are much more realistic candidates for expensive Sunny KPM calculations.
#
# Run from the repo root:
#
#   julia --project=. scripts/check_analytical_histogram_grid_convergence.jl

using Printf
using Statistics
using LinearAlgebra
using Random

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "scripts", "legacy", "YZGO_cofit_9T14T_shared_fraction_legacy.jl"))

function _get_nested(config::Dict, keys::Vector{String}, default=nothing)
    cur = config
    for key in keys
        if !(cur isa Dict) || !haskey(cur, key)
            return default
        end
        cur = cur[key]
    end
    return cur
end

_as_string_vector(x) = x === nothing ? String[] : String[string(v) for v in x]
_as_float_vector(x) = x === nothing ? Float64[] : Float64[Float64(v) for v in x]
_as_int_vector(x) = x === nothing ? Int[] : Int[round(Int, v) for v in x]

function _fit_param_dict_from_specs(specs)
    p = Dict{Symbol,Float64}()
    for s in specs
        p[s.name] = Float64(s.initial)
    end
    return p
end

_model_key_for_qtag(qtag::String) = NF.OVERLAY_QTAG_TO_MODELKEY[qtag]

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

function _midpoint_nodes(range::Tuple{Float64,Float64}, n::Integer)
    n <= 1 && return ([0.5 * (range[1] + range[2])], [1.0])
    xs = [range[1] + (i - 0.5) * (range[2] - range[1]) / n for i in 1:n]
    ws = fill(1.0 / n, n)
    return xs, ws
end

function _gaussian_grid_nodes(n::Integer; nsigma::Real=3.0)
    n <= 1 && return ([0.0], [1.0])
    zs = collect(range(-Float64(nsigma), Float64(nsigma); length=n))
    ws = exp.(-0.5 .* zs.^2)
    ws ./= sum(ws)
    return zs, ws
end

function _finite_resolution_sigmas(resolution::NF.ResolutionParams)
    if !resolution.momentum_enabled
        return (0.0, 0.0, 0.0)
    end
    resolution.momentum_mode == :rlu_diagonal || error("Grid convergence script currently supports resolution.momentum_mode = :rlu_diagonal only; got $(resolution.momentum_mode)")
    return (resolution.sigma_H_rlu, resolution.sigma_K_rlu, resolution.sigma_L_rlu)
end

function simulate_cut_1d_grid(cut::NF.CutSpec1D, model::NF.ModelParams;
                              lattice::NF.LatticeParams = NF.LatticeParams(),
                              disorder::NF.DisorderParams = NF.DisorderParams(),
                              resolution::NF.ResolutionParams = NF.ResolutionParams(),
                              rng::AbstractRNG = MersenneTwister(1234),
                              n_measured_h::Int = 3,
                              n_measured_k::Int = 3,
                              n_measured_l::Int = 1,
                              n_resolution_h::Int = 3,
                              n_resolution_k::Int = 3,
                              n_resolution_l::Int = 1,
                              resolution_nsigma::Real = 3.0,
                              n_disorder_per_q::Int = 1,
                              use_form_factor::Bool = true,
                              include_j2_formfactor::Bool = true,
                              polarization::Symbol = :transverse_c,
                              include_kfki::Bool = true,
                              intensity_scale::Real = 1.0)

    edges = collect(cut.E_range[1]:cut.dE_meV:cut.E_range[2])
    centers = NF.bin_centers(edges)
    hist = zeros(Float64, length(centers))

    Hnodes, Hw = _midpoint_nodes(cut.H_range, n_measured_h)
    Knodes, Kw = _midpoint_nodes(cut.K_range, n_measured_k)
    Lnodes, Lw = _midpoint_nodes(cut.L_range, n_measured_l)

    zH, wH = _gaussian_grid_nodes(n_resolution_h; nsigma=resolution_nsigma)
    zK, wK = _gaussian_grid_nodes(n_resolution_k; nsigma=resolution_nsigma)
    zL, wL = _gaussian_grid_nodes(n_resolution_l; nsigma=resolution_nsigma)
    sigmaH, sigmaK, sigmaL = _finite_resolution_sigmas(resolution)

    n_dis = max(1, n_disorder_per_q)
    n_q_grid = length(Hnodes) * length(Knodes) * length(Lnodes) * length(zH) * length(zK) * length(zL)
    n_events = n_q_grid * n_dis
    naccepted = 0

    for (iH, Hm) in pairs(Hnodes), (iK, Km) in pairs(Knodes), (iL, Lm) in pairs(Lnodes)
        w_meas = Hw[iH] * Kw[iK] * Lw[iL]
        for (jH, zh) in pairs(zH), (jK, zk) in pairs(zK), (jL, zl) in pairs(zL)
            Htrue = Float64(Hm) + sigmaH * Float64(zh)
            Ktrue = Float64(Km) + sigmaK * Float64(zk)
            Ltrue = Float64(Lm) + sigmaL * Float64(zl)
            w_q = w_meas * wH[jH] * wK[jK] * wL[jL]
            for _ in 1:n_dis
                gzz, J1, J2 = NF.sample_disordered_parameters(model, disorder, rng)
                Etrue = NF.dispersion_meV(Htrue, Ktrue, model; gzz=gzz, J1_meV=J1, J2_meV=J2)
                if !isfinite(Etrue) || Etrue < 0.0 || Etrue >= resolution.Ei_meV
                    continue
                end
                naccepted += 1

                qtrue = NF.rlu_to_cart(Htrue, Ktrue, Ltrue, lattice)
                Qmag = norm(qtrue)
                ff = use_form_factor ? NF.yb3_form_factor(Qmag; include_j2=include_j2_formfactor) : 1.0
                pol = NF.polarization_factor(qtrue; mode=polarization)
                kin = include_kfki ? NF.kf_over_ki(Etrue, resolution.Ei_meV) : 1.0
                weight = Float64(intensity_scale) * w_q / n_dis * ff^2 * pol * kin

                if resolution.energy_enabled
                    sigmaE = NF.energy_resolution_sigma_meV(Etrue, resolution.Ei_meV)
                    NF.deposit_gaussian_1d!(hist, edges, Etrue, sigmaE, weight)
                else
                    NF.deposit_delta_1d!(hist, edges, Etrue, weight)
                end
            end
        end
    end

    return (
        name = cut.name,
        cut = cut,
        model = model,
        lattice = lattice,
        disorder = disorder,
        resolution = resolution,
        n_samples = n_events,
        n_q_grid = n_q_grid,
        n_disorder_per_q = n_dis,
        naccepted = naccepted,
        acceptance_fraction = n_events > 0 ? naccepted / n_events : NaN,
        E_edges_meV = edges,
        E_centers_meV = centers,
        intensity = hist,
        sampling_mode = :deterministic_measured_Q_resolution_grid,
    )
end

function simulate_all_cuts_1d_grid(model::NF.ModelParams;
                                   cuts::Vector{NF.CutSpec1D} = NF.default_cuts_1d(),
                                   kwargs...)
    out = Dict{String, Any}()
    for cut in cuts
        out[cut.name] = simulate_cut_1d_grid(cut, model; kwargs...)
    end
    return out
end

function _legacy_sigma(p::Dict{Symbol,Float64}, nm::Symbol, fallback::Symbol)
    return haskey(p, nm) ? p[nm] : get(p, fallback, 0.0)
end

function simulate_two_kernel_model_fields_grid(p::Dict{Symbol,Float64};
                                               fields_T::AbstractVector{<:Real}=[9.0],
                                               cuts::Vector{NF.CutSpec1D}=NF.default_cuts_1d(),
                                               lattice::NF.LatticeParams=NF.demo_defaults().lattice,
                                               resolution::NF.ResolutionParams=NF.demo_defaults().resolution,
                                               seed::Int=2026,
                                               S::Real=0.5,
                                               J_units::Symbol=:fractional,
                                               correlate_J1_J2::Bool=false,
                                               use_form_factor::Bool=true,
                                               include_j2_formfactor::Bool=true,
                                               include_kfki::Bool=true,
                                               polarization::Symbol=:transverse_c,
                                               n_measured_h::Int=3,
                                               n_measured_k::Int=3,
                                               n_measured_l::Int=1,
                                               n_resolution_h::Int=3,
                                               n_resolution_k::Int=3,
                                               n_resolution_l::Int=1,
                                               resolution_nsigma::Real=3.0,
                                               n_disorder_per_q::Int=1)
    disorder_disp = NF.DisorderParams(
        enabled = true,
        sigma_gzz = max(p[:sigma_gzz], 0.0),
        sigma_J1 = max(_legacy_sigma(p, :sigma_J1, :sigma_J), 0.0),
        sigma_J2 = max(_legacy_sigma(p, :sigma_J2, :sigma_J), 0.0),
        J_units = J_units,
        correlate_J1_J2 = correlate_J1_J2,
    )
    disorder_flat = NF.DisorderParams(
        enabled = true,
        sigma_gzz = max(p[:sigma_gzz2], 0.0),
        sigma_J1 = 0.0,
        sigma_J2 = 0.0,
        J_units = :absolute_meV,
        correlate_J1_J2 = false,
    )
    gperp_scales = NF.neutron_gperp_intensity_scales(p)

    out = Dict{Float64, Any}()
    for B in fields_T
        Bf = Float64(B)
        model_disp = NF.ModelParams(B_T=Bf, gzz=p[:gzz], J1_meV=p[:J1_meV], J2_meV=p[:J2_meV], S=Float64(S))
        model_flat = NF.ModelParams(B_T=Bf, gzz=p[:gzz2], J1_meV=0.0, J2_meV=0.0, S=Float64(S))

        field_seed = seed + round(Int, 1000 * Bf)
        rng_disp = MersenneTwister(field_seed)
        rng_flat = MersenneTwister(field_seed + 7919)

        common_kwargs = (
            cuts = cuts,
            lattice = lattice,
            resolution = resolution,
            n_measured_h = n_measured_h,
            n_measured_k = n_measured_k,
            n_measured_l = n_measured_l,
            n_resolution_h = n_resolution_h,
            n_resolution_k = n_resolution_k,
            n_resolution_l = n_resolution_l,
            resolution_nsigma = resolution_nsigma,
            n_disorder_per_q = n_disorder_per_q,
            use_form_factor = use_form_factor,
            include_j2_formfactor = include_j2_formfactor,
            include_kfki = include_kfki,
            polarization = polarization,
        )

        dispersive = simulate_all_cuts_1d_grid(
            model_disp;
            disorder = disorder_disp,
            rng = rng_disp,
            intensity_scale = gperp_scales.dispersive,
            common_kwargs...,
        )
        nondispersive = simulate_all_cuts_1d_grid(
            model_flat;
            disorder = disorder_flat,
            rng = rng_flat,
            intensity_scale = gperp_scales.nondispersive,
            common_kwargs...,
        )
        out[Bf] = (dispersive=dispersive, nondispersive=nondispersive)
    end
    return out
end

function evaluate_neutron_fit_grid(p::Dict{Symbol,Float64}, data_scans;
                                   fields_T::AbstractVector{<:Real}=[9.0],
                                   qtags=["0_1_0", "0p33_0p33_0", "0p5_0_0"],
                                   fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
                                   seed::Int=2026,
                                   use_errors::Bool=true,
                                   error_floor::Real=0.0,
                                   lattice::NF.LatticeParams=NF.demo_defaults().lattice,
                                   resolution::NF.ResolutionParams=NF.demo_defaults().resolution,
                                   cuts::Vector{NF.CutSpec1D}=NF.default_cuts_1d(),
                                   S::Real=0.5,
                                   J_units::Symbol=:fractional,
                                   correlate_J1_J2::Bool=false,
                                   use_form_factor::Bool=true,
                                   include_j2_formfactor::Bool=true,
                                   include_kfki::Bool=true,
                                   polarization::Symbol=:transverse_c,
                                   nfree::Int=1,
                                   n_measured_h::Int=3,
                                   n_measured_k::Int=3,
                                   n_measured_l::Int=1,
                                   n_resolution_h::Int=3,
                                   n_resolution_k::Int=3,
                                   n_resolution_l::Int=1,
                                   resolution_nsigma::Real=3.0,
                                   n_disorder_per_q::Int=1)
    pn = neutron_param_dict_from_cofit(p)
    r2n = NF.second_kernel_relative_intensity(pn)

    model = simulate_two_kernel_model_fields_grid(pn;
        fields_T = fields_T,
        cuts = cuts,
        lattice = lattice,
        resolution = resolution,
        seed = seed,
        S = S,
        J_units = J_units,
        correlate_J1_J2 = correlate_J1_J2,
        use_form_factor = use_form_factor,
        include_j2_formfactor = include_j2_formfactor,
        include_kfki = include_kfki,
        polarization = polarization,
        n_measured_h = n_measured_h,
        n_measured_k = n_measured_k,
        n_measured_l = n_measured_l,
        n_resolution_h = n_resolution_h,
        n_resolution_k = n_resolution_k,
        n_resolution_l = n_resolution_l,
        resolution_nsigma = resolution_nsigma,
        n_disorder_per_q = n_disorder_per_q,
    )

    bundle = NF.build_fit_vector_bundle(data_scans, model;
        fields_T = fields_T,
        qtags = qtags,
        fit_windows_by_q = fit_windows_by_q,
        use_errors = use_errors,
        error_floor = error_floor,
    )

    if length(bundle.data_intensity) < 6
        return (; redchi2=Inf, scale=NaN, r2=r2n, model=model, bundle=bundle)
    end

    scale = if haskey(p, :log10_neutron_scale)
        10.0 ^ p[:log10_neutron_scale]
    else
        NF.best_global_scale(bundle, r2n)
    end
    redchi2 = NF.chisq_for_bundle(bundle, scale, r2n; nfree=nfree, reduced=true)
    return (; redchi2=redchi2, scale=scale, r2=r2n, model=model, bundle=bundle)
end

function _write_grid_inventory(path::AbstractString, rows)
    header = [
        "label", "n_measured_h", "n_measured_k", "n_measured_l",
        "n_resolution_h", "n_resolution_k", "n_resolution_l", "resolution_nsigma",
        "n_disorder_per_q", "q_grid_points_per_cut", "component_events_per_cut", "total_component_events_6cuts",
    ]
    _write_csv(path, header, rows)
end

function _write_metrics(path::AbstractString; rows)
    header = [
        "label", "reference_label", "field_T", "qtag", "component",
        "rmse", "relative_rmse", "max_abs_diff", "integrated", "reference_integrated",
        "relative_integrated_error", "peak_E_meV", "reference_peak_E_meV", "peak_E_shift_meV",
        "peak_I", "reference_peak_I", "relative_peak_height_error",
        "linewidth_sigma_meV", "reference_linewidth_sigma_meV", "relative_linewidth_error",
        "q_grid_points_per_cut", "component_events_per_cut",
    ]
    _write_csv(path, header, rows)
end

function _write_model_long_csv(path::AbstractString, results_by_label::Dict{String,Any}, grid_meta::Dict{String,Any}, fields_T, qtags)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "label,field_T,qtag,energy_meV,I_total_scaled,I_dispersive_scaled,I_flat_scaled,I_dispersive_unscaled,I_flat_unscaled,neutron_global_scale,r2_shared,q_grid_points_per_cut,component_events_per_cut")
        for label in sort(collect(keys(results_by_label)))
            nres = results_by_label[label]
            meta = grid_meta[label]
            for B in fields_T, qtag in qtags
                arr = _component_arrays(nres, B, qtag)
                for i in eachindex(arr.E)
                    println(io, join((
                        label,
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
                        meta.q_grid_points_per_cut,
                        meta.component_events_per_cut,
                    ), ","))
                end
            end
        end
    end
    return path
end

function _make_grid_convergence_plots(figdir::AbstractString, results_by_label::Dict{String,Any}, metrics_rows, fields_T, qtags, overlay_labels)
    try
        @eval import CairoMakie as CM
    catch err
        @warn "Could not load CairoMakie; skipping grid convergence plots" exception=(err, catch_backtrace())
        return String[]
    end
    mkpath(figdir)
    files = String[]
    labels = [lab for lab in overlay_labels if haskey(results_by_label, lab)]
    isempty(labels) && return String[]

    fig = CM.Figure(size=(1500, 800))
    for (ir, B) in enumerate(fields_T)
        for (ic, qtag) in enumerate(qtags)
            ax = CM.Axis(fig[ir, ic];
                title=@sprintf("%s, %.0f T", qtag, Float64(B)),
                xlabel="Energy transfer (meV)",
                ylabel=ic == 1 ? "Analytical model intensity" : "",
            )
            for lab in labels
                arr = _component_arrays(results_by_label[lab], B, qtag)
                lw = occursin("reference", lowercase(lab)) || occursin("mc", lowercase(lab)) ? 3.0 : 1.5
                CM.lines!(ax, arr.E, arr.total; label=lab, linewidth=lw)
            end
            if ir == 1 && ic == length(qtags)
                CM.axislegend(ax; position=:rt, framevisible=false)
            end
        end
    end
    path_overlay = joinpath(figdir, "analytical_grid_convergence_1d_overlays.png")
    CM.save(path_overlay, fig)
    push!(files, path_overlay)

    fig2 = CM.Figure(size=(1100, 700))
    ax2 = CM.Axis(fig2[1, 1];
        title="Analytical deterministic-grid convergence vs MC reference",
        xlabel="component-events per cut",
        ylabel="relative RMSE vs reference",
        xscale=log10,
    )
    # Metric tuple layout: label, reference_label, field, qtag, component, ..., relative_rmse=row[7], component_events=row[22]
    for B in fields_T, qtag in qtags
        xs = Float64[]
        ys = Float64[]
        for row in metrics_rows
            if Float64(row[3]) == Float64(B) && row[4] == qtag && row[5] == "total"
                push!(xs, Float64(row[22]))
                push!(ys, Float64(row[7]))
            end
        end
        if !isempty(xs)
            order = sortperm(xs)
            CM.lines!(ax2, xs[order], ys[order]; label=@sprintf("%s %.0fT", qtag, Float64(B)), linewidth=2)
            CM.scatter!(ax2, xs[order], ys[order]; markersize=9)
        end
    end
    CM.axislegend(ax2; position=:rt, framevisible=false)
    path_metrics = joinpath(figdir, "analytical_grid_convergence_relative_rmse.png")
    CM.save(path_metrics, fig2)
    push!(files, path_metrics)
    return files
end

function _parse_grid_configs(controls::Dict)
    raw = _get_nested(controls, ["grid", "configs"], nothing)
    if raw === nothing
        return Any[]
    end
    return raw
end

function main()
    controls_path = joinpath(REPO_ROOT, "configs", "analytical_histogram_grid_convergence_controls.toml")
    cofit_controls_path = joinpath(REPO_ROOT, "configs", "cofit_controls.toml")
    best_fit_path = joinpath(REPO_ROOT, "configs", "best_fit_parameters.toml")

    controls = load_toml_config(controls_path)
    cofit_controls = load_cofit_controls(cofit_controls_path)
    best_fit_config = load_best_fit_parameters(best_fit_path)

    initial_guess_kwargs = cofit_initial_guess_kwargs(best_fit_config)
    specs = cofit_default_param_specs(; initial_guess_kwargs...)
    pbest = _fit_param_dict_from_specs(specs)

    table_dir = joinpath(REPO_ROOT, string(_get_nested(controls, ["output", "table_subdir"], "results/feature_tables/analytical_histogram_grid_convergence")))
    fig_dir = joinpath(REPO_ROOT, string(_get_nested(controls, ["output", "figure_subdir"], "results/figures/analytical_histogram_grid_convergence")))
    mkpath(table_dir)
    mkpath(fig_dir)

    fields_T = _as_float_vector(cofit_controls["data"]["fields_T"])
    qtags = _as_string_vector(cofit_controls["data"]["qtags"])
    base_dir = joinpath(REPO_ROOT, "data", "neutron", "CNCS_1d_scans")
    data_mode = toml_symbol(cofit_controls["data"]["data_mode"])

    model_controls = controls["model"]
    S = Float64(get(model_controls, "S", 0.5))
    J_units = toml_symbol(get(model_controls, "J_units", "fractional"))
    correlate_J1_J2 = Bool(get(model_controls, "correlate_J1_J2", false))
    use_form_factor = Bool(get(model_controls, "use_form_factor", true))
    include_j2_formfactor = Bool(get(model_controls, "include_j2_formfactor", true))
    include_kfki = Bool(get(model_controls, "include_kfki", true))
    polarization = toml_symbol(get(model_controls, "polarization", "transverse_c"))

    seed = round(Int, _get_nested(controls, ["sampling", "seed"], 2026))
    reference_n = round(Int, _get_nested(controls, ["reference", "mc_n_samples_per_cut"], 500000))
    reference_label = string(_get_nested(controls, ["reference", "label"], @sprintf("MC_%d", reference_n)))
    grid_configs = _parse_grid_configs(controls)
    isempty(grid_configs) && error("No [[grid.configs]] entries found in $controls_path")

    println("Analytical deterministic-grid histogram convergence study")
    println("--------------------------------------------------------")
    println("controls:        ", controls_path)
    println("best parameters: ", best_fit_path)
    println("data mode:       ", data_mode)
    println("reference:       ", reference_label, " with N=", reference_n, " MC samples/cut/component")
    println("fields:          ", fields_T)
    println("qtags:           ", qtags)
    println("grid configs:    ", length(grid_configs))
    println("output tables:   ", table_dir)
    println("output figures:  ", fig_dir)
    println()

    data_scans, _, _ = load_neutron_fit_data_1d_filtered(;
        base_dir=base_dir,
        data_mode=data_mode,
        Ei_meV=cofit_controls["data"]["neutron_fit_Ei_meV"],
        temperature_K=cofit_controls["data"]["neutron_fit_temperature_K"],
    )

    println(@sprintf("Running MC reference with n_samples_per_cut = %d", reference_n))
    ref = evaluate_neutron_fit(pbest, data_scans;
        fields_T=fields_T,
        qtags=qtags,
        fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
        n_samples_per_cut=reference_n,
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
    println(@sprintf("  reference redchi2 = %.8g, scale = %.8g, r2 = %.8g", ref.redchi2, ref.scale, ref.r2))

    results_by_label = Dict{String,Any}(reference_label => ref)
    grid_meta = Dict{String,Any}(reference_label => (; q_grid_points_per_cut=reference_n, component_events_per_cut=reference_n))
    inventory_rows = Tuple[]

    for cfg in grid_configs
        label = string(get(cfg, "label", "grid"))
        n_measured_h = round(Int, get(cfg, "n_measured_h", 3))
        n_measured_k = round(Int, get(cfg, "n_measured_k", 3))
        n_measured_l = round(Int, get(cfg, "n_measured_l", 1))
        n_resolution_h = round(Int, get(cfg, "n_resolution_h", 3))
        n_resolution_k = round(Int, get(cfg, "n_resolution_k", 3))
        n_resolution_l = round(Int, get(cfg, "n_resolution_l", 1))
        resolution_nsigma = Float64(get(cfg, "resolution_nsigma", 3.0))
        n_disorder_per_q = round(Int, get(cfg, "n_disorder_per_q", 1))
        q_grid_points_per_cut = n_measured_h * n_measured_k * n_measured_l * n_resolution_h * n_resolution_k * n_resolution_l
        component_events_per_cut = q_grid_points_per_cut * max(1, n_disorder_per_q)

        println()
        println(@sprintf("Running grid config %-20s: measured=(%d,%d,%d), resolution=(%d,%d,%d), disorder/q=%d, q_grid=%d, component_events/cut=%d",
            label, n_measured_h, n_measured_k, n_measured_l, n_resolution_h, n_resolution_k, n_resolution_l,
            n_disorder_per_q, q_grid_points_per_cut, component_events_per_cut))

        nres = evaluate_neutron_fit_grid(pbest, data_scans;
            fields_T=fields_T,
            qtags=qtags,
            fit_windows_by_q=NF.YZGO_FIT_WINDOWS_1D,
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
            n_measured_h=n_measured_h,
            n_measured_k=n_measured_k,
            n_measured_l=n_measured_l,
            n_resolution_h=n_resolution_h,
            n_resolution_k=n_resolution_k,
            n_resolution_l=n_resolution_l,
            resolution_nsigma=resolution_nsigma,
            n_disorder_per_q=n_disorder_per_q,
        )
        results_by_label[label] = nres
        grid_meta[label] = (; q_grid_points_per_cut=q_grid_points_per_cut, component_events_per_cut=component_events_per_cut)
        push!(inventory_rows, (
            label, n_measured_h, n_measured_k, n_measured_l,
            n_resolution_h, n_resolution_k, n_resolution_l, resolution_nsigma,
            n_disorder_per_q, q_grid_points_per_cut, component_events_per_cut, 6 * component_events_per_cut,
        ))
        println(@sprintf("  redchi2 = %.8g, scale = %.8g, r2 = %.8g", nres.redchi2, nres.scale, nres.r2))
    end

    components = _as_string_vector(_get_nested(controls, ["comparison", "components"], ["total", "dispersive", "flat"]))
    metrics_rows = Tuple[]
    for (label, nres) in results_by_label
        label == reference_label && continue
        meta = grid_meta[label]
        for B in fields_T, qtag in qtags
            arr = _component_arrays(nres, B, qtag)
            arr_ref = _component_arrays(ref, B, qtag)
            for comp in components
                y = getproperty(arr, Symbol(comp))
                yref = getproperty(arr_ref, Symbol(comp))
                m = _peak_metrics(arr.E, y)
                mref = _peak_metrics(arr_ref.E, yref)
                rmse = _rmse(y, yref)
                relrmse = _relative_rmse(y, yref)
                maxabs = maximum(abs.(Float64.(y) .- Float64.(yref)))
                relint = (m.integrated - mref.integrated) / max(abs(mref.integrated), eps(Float64))
                relpeak = (m.peak_I - mref.peak_I) / max(abs(mref.peak_I), eps(Float64))
                rellw = (m.linewidth_sigma - mref.linewidth_sigma) / max(abs(mref.linewidth_sigma), eps(Float64))
                push!(metrics_rows, (
                    label, reference_label, Float64(B), qtag, comp,
                    rmse, relrmse, maxabs,
                    m.integrated, mref.integrated, relint,
                    m.peak_E, mref.peak_E, m.peak_E - mref.peak_E,
                    m.peak_I, mref.peak_I, relpeak,
                    m.linewidth_sigma, mref.linewidth_sigma, rellw,
                    meta.q_grid_points_per_cut, meta.component_events_per_cut,
                ))
            end
        end
    end

    inventory_path = joinpath(table_dir, "analytical_grid_convergence_inventory.csv")
    metrics_path = joinpath(table_dir, "analytical_grid_convergence_metrics_vs_mc_reference.csv")
    long_path = joinpath(table_dir, "analytical_grid_convergence_models_long.csv")
    _write_grid_inventory(inventory_path, inventory_rows)
    _write_metrics(metrics_path; rows=metrics_rows)
    _write_model_long_csv(long_path, results_by_label, grid_meta, fields_T, qtags)

    println()
    println("Wrote grid convergence tables:")
    println("  ", inventory_path)
    println("  ", metrics_path)
    println("  ", long_path)

    if Bool(_get_nested(controls, ["plotting", "make_plots"], true))
        overlay_labels = _as_string_vector(_get_nested(controls, ["plotting", "overlay_labels"], [reference_label, "grid_3x3_meas_3x3_res", "grid_5x5_meas_5x5_res"]))
        files = _make_grid_convergence_plots(fig_dir, results_by_label, metrics_rows, fields_T, qtags, overlay_labels)
        if !isempty(files)
            println("Wrote grid convergence figures:")
            foreach(f -> println("  ", f), files)
        end
    end

    println()
    println("Grid convergence study complete.")
    return (; results_by_label=results_by_label, metrics_path=metrics_path, long_path=long_path, inventory_path=inventory_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
