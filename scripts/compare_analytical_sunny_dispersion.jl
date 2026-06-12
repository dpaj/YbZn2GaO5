# scripts/compare_analytical_sunny_dispersion.jl
#
# Clean diagnostic comparison between the analytical field-polarized
# YbZn2GaO5 triangular-lattice dispersion and Sunny.jl LSWT.
#
# This script deliberately avoids disorder, KPM, flat/nondispersive components,
# experimental histogramming, and fitted neutron scale factors. It is only a
# clean geometry / exchange-shell / dispersion check for the dispersive model.
#
# Run from the repo root:
#   julia --project=. scripts/compare_analytical_sunny_dispersion.jl
#
# Convention used here:
#   The effective triangular lattice is built as a P1 one-site lattice with
#   gamma = 120 degrees. Because P1 has no rotational symmetry, all three
#   positive representatives of each exchange shell are set explicitly:
#       J1: [1,0,0], [0,1,0], [1,1,0]
#       J2: [1,-1,0], [2,1,0], [1,2,0]
#   These are exactly the bond shells corresponding to the analytical form
#   factors
#       Δ1 = 6 - 2[cos(2πH) + cos(2πK) + cos(2π(H+K))]
#       Δ2 = 6 - 2[cos(2π(H-K)) + cos(2π(2H+K)) + cos(2π(H+2K))]
#

using LinearAlgebra
using Printf
using StaticArrays
using DelimitedFiles
using Sunny
using CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "parameters.jl"))

const MU_B_MEV_PER_T = 0.05788381806

# -----------------------------------------------------------------------------
# Small local helpers
# -----------------------------------------------------------------------------

function get_nested(d::Dict, keys::Vector{String}, default)
    cur = d
    for (i, k) in enumerate(keys)
        if !(cur isa Dict) || !haskey(cur, k)
            return default
        end
        if i == length(keys)
            return cur[k]
        end
        cur = cur[k]
    end
    return default
end

function repo_path(rel_or_abs::AbstractString)
    s = String(rel_or_abs)
    if isabspath(s)
        return normpath(s)
    end
    parts = filter(!isempty, split(s, r"[\\/]+"))
    return normpath(joinpath(REPO_ROOT, parts...))
end

function ensure_dir(path::AbstractString)
    mkpath(path)
    return path
end

function csv_cell(x)
    if x isa AbstractString
        return replace(x, "," => ";")
    elseif x isa Symbol
        return String(x)
    elseif x isa Real
        y = Float64(x)
        if isnan(y)
            return "NaN"
        elseif isinf(y)
            return y > 0 ? "Inf" : "-Inf"
        else
            return string(x)
        end
    else
        return string(x)
    end
end

function write_namedtuple_csv(path::AbstractString, rows::Vector{<:NamedTuple})
    open(path, "w") do io
        if isempty(rows)
            return
        end
        names = propertynames(first(rows))
        println(io, join(String.(names), ","))
        for r in rows
            println(io, join((csv_cell(getproperty(r, n)) for n in names), ","))
        end
    end
end

# -----------------------------------------------------------------------------
# Analytical polarized-state dispersion and exchange form factors
# -----------------------------------------------------------------------------

const J1_BONDS_ANALYTICAL_120 = NTuple{3,Int}[(1,0,0), (0,1,0), (1,1,0)]
const J2_BONDS_ANALYTICAL_120 = NTuple{3,Int}[(1,-1,0), (2,1,0), (1,2,0)]

function delta_from_bonds(H::Real, K::Real, L::Real, bonds::Vector{NTuple{3,Int}})
    Hs = Float64(H)
    Ks = Float64(K)
    Ls = Float64(L)
    return sum(2.0 * (1.0 - cos(2.0*pi*(b[1]*Hs + b[2]*Ks + b[3]*Ls))) for b in bonds)
end

function delta1_triangular(H::Real, K::Real)
    return 6.0 - 2.0 * (
        cos(2.0*pi*H) +
        cos(2.0*pi*K) +
        cos(2.0*pi*(H + K))
    )
end

function delta2_triangular(H::Real, K::Real)
    return 6.0 - 2.0 * (
        cos(2.0*pi*(H - K)) +
        cos(2.0*pi*(2.0*H + K)) +
        cos(2.0*pi*(H + 2.0*K))
    )
end

function analytical_dispersion_meV(H::Real, K::Real, params; field_T::Real, S::Real=0.5)
    h = params.gzz * MU_B_MEV_PER_T * Float64(field_T)
    return h - Float64(S) * (
        params.J1_meV * delta1_triangular(H, K) +
        params.J2_meV * delta2_triangular(H, K)
    )
end

# -----------------------------------------------------------------------------
# Sunny clean LSWT builder
# -----------------------------------------------------------------------------

function effective_triangle_crystal(controls::Dict)
    a = Float64(get_nested(controls, ["common", "lattice_a_angstrom"], 3.4))
    c = Float64(get_nested(controls, ["common", "lattice_c_angstrom"], 10.0))

    # P1 is intentional. It avoids hidden space-group assumptions and makes the
    # explicit bond list below the single source of truth for the shell geometry.
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    return Crystal(latvecs, [[0.0, 0.0, 0.0]], 1; types=["Yb"])
end

function set_g_tensor_all_sites!(sys, gzz::Real; gxy::Real=1.0)
    for site in eachsite(sys)
        Gm = zeros(Float64, 3, 3)
        Gm[1,1] = Float64(gxy)
        Gm[2,2] = Float64(gxy)
        Gm[3,3] = Float64(gzz)
        sys.gs[site] = SMatrix{3,3,Float64,9}(Gm)
    end
    return sys
end

function field_direction(controls::Dict)
    u = Float64.(get_nested(controls, ["common", "field_direction"], [0.0, 0.0, 1.0]))
    n = norm(u)
    n > 0 || error("[common].field_direction must be nonzero")
    return u ./ n
end

function add_explicit_exchange_shell!(sys, J::Real, bonds::Vector{NTuple{3,Int}})
    for b in bonds
        set_exchange!(sys, Float64(J), Bond(1, 1, [b[1], b[2], b[3]]))
    end
    return sys
end

function build_clean_sunny_dispersive_system(params, controls::Dict; field_T::Real)
    units = Units(:meV, :angstrom)
    cryst = effective_triangle_crystal(controls)
    S = Float64(get_nested(controls, ["common", "spin_S"], 0.5))
    seed = Int(get_nested(controls, ["common", "seed"], 20260611))
    dims_raw = get_nested(controls, ["dispersion_check", "dims"], [1, 1, 1])
    dims = (Int(dims_raw[1]), Int(dims_raw[2]), Int(dims_raw[3]))
    gxy = Float64(get_nested(controls, ["common", "sunny_transverse_gxy"], 1.0))

    sys = System(cryst, [1 => Moment(s=S, g=Float64(params.gzz))], :dipole; dims, seed)
    set_g_tensor_all_sites!(sys, params.gzz; gxy)

    # Explicit P1 bond shell convention matching the analytical form factors.
    add_explicit_exchange_shell!(sys, params.J1_meV, J1_BONDS_ANALYTICAL_120)
    add_explicit_exchange_shell!(sys, params.J2_meV, J2_BONDS_ANALYTICAL_120)

    u = field_direction(controls)
    set_field!(sys, collect(u .* (Float64(field_T) * units.T)))

    # At 9/14 T this should be the field-polarized state. Use the standard
    # randomize/minimize route so this script stays close to the Sunny examples.
    randomize_spins!(sys)
    maxiters = Int(get_nested(controls, ["dispersion_check", "maxiters"], 1000))
    minimize_energy!(sys; maxiters)

    return (; sys, crystal=cryst, units, gxy, dims)
end

# -----------------------------------------------------------------------------
# Path construction and Sunny result extraction
# -----------------------------------------------------------------------------

function default_path_vertices(controls::Dict)
    verts = get_nested(controls, ["dispersion_check", "path_vertices_exp"], nothing)
    if verts === nothing
        return [[0.0, 0.0, 0.0], [1/3, 1/3, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]]
    end
    return [Float64.(v) for v in verts]
end

function vertex_labels(controls::Dict, nverts::Int)
    labels = get_nested(controls, ["dispersion_check", "path_labels"], nothing)
    if labels === nothing
        defaults = ["Γ", "K", "M", "Γ"]
        return defaults[1:min(nverts, length(defaults))]
    end
    return String.(labels)
end

function cumulative_distance(qs::AbstractVector)
    x = zeros(Float64, length(qs))
    for i in 2:length(qs)
        dq = Float64.(qs[i]) .- Float64.(qs[i-1])
        x[i] = x[i-1] + norm(dq)
    end
    return x
end

function extract_dispersion_and_intensity(res, nq::Int)
    hasproperty(res, :disp) || error("Sunny intensities_bands result has no .disp field. propertynames(res)=$(propertynames(res))")
    disp = Float64.(getproperty(res, :disp))
    ndims(disp) == 2 || error("Expected res.disp to be a 2D array; got size $(size(disp))")

    if size(disp, 2) == nq
        energies = disp
    elseif size(disp, 1) == nq
        energies = permutedims(disp)
    else
        error("Could not align Sunny res.disp size $(size(disp)) with nq=$nq")
    end

    data = fill(NaN, size(energies))
    if hasproperty(res, :data)
        raw = Float64.(getproperty(res, :data))
        if ndims(raw) == 2 && size(raw) == size(energies)
            data .= raw
        elseif ndims(raw) == 2 && size(permutedims(raw)) == size(energies)
            data .= permutedims(raw)
        end
    end
    return (; energies, intensities=data)
end

function compute_sunny_dispersion(params, controls::Dict; field_T::Real)
    built = build_clean_sunny_dispersive_system(params, controls; field_T)
    qverts = default_path_vertices(controls)
    npts = Int(get_nested(controls, ["dispersion_check", "n_path_points"], 301))

    path = q_space_path(built.crystal, qverts, npts)
    swt = SpinWaveTheory(built.sys; measure=ssf_perp(built.sys))
    res = intensities_bands(swt, path)
    bands = extract_dispersion_and_intensity(res, length(path.qs))

    qs = [Float64.(q) for q in path.qs]
    x = cumulative_distance(qs)
    S = Float64(get_nested(controls, ["common", "spin_S"], 0.5))
    Ean = [analytical_dispersion_meV(q[1], q[2], params; field_T, S) for q in qs]

    tick_qs = [Float64.(q) for q in qverts]
    tick_x = cumulative_distance(tick_qs)
    labels = vertex_labels(controls, length(tick_x))

    return (; field_T=Float64(field_T), path, qs, x,
        E_analytical=Ean, E_sunny=bands.energies, I_sunny=bands.intensities,
        tick_x, tick_labels=labels, gxy=built.gxy, dims=built.dims)
end

function nearest_path_index(qs, target)
    d = [norm(q .- target) for q in qs]
    return argmin(d)
end

# -----------------------------------------------------------------------------
# Geometry sanity table
# -----------------------------------------------------------------------------

function make_geometry_rows()
    points = [
        (label="Γ", H=0.0, K=0.0, L=0.0),
        (label="K_1over3_1over3", H=1/3, K=1/3, L=0.0),
        (label="M_0p5_0_0", H=0.5, K=0.0, L=0.0),
        (label="Γ_equiv_1_0_0", H=1.0, K=0.0, L=0.0),
    ]
    rows = NamedTuple[]
    for p in points
        d1a = delta1_triangular(p.H, p.K)
        d2a = delta2_triangular(p.H, p.K)
        d1b = delta_from_bonds(p.H, p.K, p.L, J1_BONDS_ANALYTICAL_120)
        d2b = delta_from_bonds(p.H, p.K, p.L, J2_BONDS_ANALYTICAL_120)
        push!(rows, (; p.label, p.H, p.K, p.L,
            delta1_analytical=d1a, delta1_bond_shell=d1b, delta1_difference=d1b-d1a,
            delta2_analytical=d2a, delta2_bond_shell=d2b, delta2_difference=d2b-d2a))
    end
    return rows
end

function print_geometry_table(rows)
    println("Exchange-shell form-factor sanity check")
    println("---------------------------------------")
    @printf("%-18s %12s %12s %12s %12s\n", "point", "Δ1 analytical", "Δ1 bonds", "Δ2 analytical", "Δ2 bonds")
    for r in rows
        @printf("%-18s %12.6g %12.6g %12.6g %12.6g\n",
            r.label, r.delta1_analytical, r.delta1_bond_shell,
            r.delta2_analytical, r.delta2_bond_shell)
    end
    println()
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

function main()
    controls_path = joinpath(REPO_ROOT, "configs", "sunny_validation_controls.toml")
    params_path = joinpath(REPO_ROOT, "configs", "best_fit_parameters.toml")
    controls = isfile(controls_path) ? load_toml_config(controls_path) : Dict{String,Any}()
    params = load_canonical_model_parameters(params_path)

    println("Canonical model parameters used for analytical and Sunny clean LSWT")
    println("------------------------------------------------------------------")
    print_canonical_model_parameters(params)

    fields = Float64.(get_nested(controls, ["dispersion_check", "fields_T"], get_nested(controls, ["common", "fields_T"], [9.0, 14.0])))
    out_fig_dir = ensure_dir(repo_path(String(get_nested(controls, ["paths", "figure_subdir"], "results/figures/sunny_validation"))))
    out_tab_dir = ensure_dir(repo_path(String(get_nested(controls, ["paths", "table_subdir"], "results/feature_tables/sunny_validation"))))

    println("Sunny dispersion check settings")
    println("-------------------------------")
    println("fields_T      = ", fields)
    println("n_path_points = ", Int(get_nested(controls, ["dispersion_check", "n_path_points"], 301)))
    println("path vertices = ", default_path_vertices(controls))
    println("J1 bonds      = ", J1_BONDS_ANALYTICAL_120)
    println("J2 bonds      = ", J2_BONDS_ANALYTICAL_120)
    println("output dirs   = ", out_fig_dir, " ; ", out_tab_dir)
    println()

    geom_rows = make_geometry_rows()
    print_geometry_table(geom_rows)
    geom_csv = joinpath(out_tab_dir, "analytical_vs_sunny_exchange_geometry_check.csv")
    write_namedtuple_csv(geom_csv, geom_rows)

    results = Dict{Float64,Any}()
    long_rows = NamedTuple[]
    summary_rows = NamedTuple[]

    special_points = [
        ("Γ", [0.0, 0.0, 0.0]),
        ("K_1over3_1over3", [1/3, 1/3, 0.0]),
        ("M_0p5_0_0", [0.5, 0.0, 0.0]),
        ("Γ_equiv_1_0_0", [1.0, 0.0, 0.0]),
    ]

    for B in fields
        @info "Computing clean Sunny LSWT dispersion" field_T=B
        r = compute_sunny_dispersion(params, controls; field_T=B)
        results[B] = r

        for i in eachindex(r.x)
            H, K, L = r.qs[i]
            push!(long_rows, (; field_T=B, source="analytical", band_index=1,
                path_index=i, x_path=r.x[i], H, K, L, energy_meV=r.E_analytical[i], intensity=NaN))
            for b in 1:size(r.E_sunny, 1)
                push!(long_rows, (; field_T=B, source="sunny", band_index=b,
                    path_index=i, x_path=r.x[i], H, K, L, energy_meV=r.E_sunny[b,i], intensity=r.I_sunny[b,i]))
            end
        end

        for (label, qpt) in special_points
            j = nearest_path_index(r.qs, qpt)
            Ean = analytical_dispersion_meV(qpt[1], qpt[2], params; field_T=B,
                S=Float64(get_nested(controls, ["common", "spin_S"], 0.5)))
            for b in 1:size(r.E_sunny, 1)
                Es = r.E_sunny[b,j]
                push!(summary_rows, (; field_T=B, label, H=qpt[1], K=qpt[2], L=qpt[3],
                    E_analytical_meV=Ean, band_index=b,
                    E_sunny_meV=Es, difference_meV=Es-Ean))
            end
        end
    end

    long_csv = joinpath(out_tab_dir, "analytical_vs_sunny_clean_dispersion_long.csv")
    summary_csv = joinpath(out_tab_dir, "analytical_vs_sunny_clean_dispersion_summary.csv")
    write_namedtuple_csv(long_csv, long_rows)
    write_namedtuple_csv(summary_csv, summary_rows)

    fig = Figure(size=(1100, 420*length(fields)))
    for (iB, B) in enumerate(fields)
        r = results[B]
        ax = Axis(fig[iB, 1], xlabel="Path coordinate", ylabel="Energy (meV)",
            title=@sprintf("Clean dispersion comparison, B = %.3g T", B),
            xticks=(r.tick_x, r.tick_labels))
        lines!(ax, r.x, r.E_analytical; linewidth=4, label="analytical clean")

        for b in 1:size(r.E_sunny, 1)
            lab = @sprintf("Sunny LSWT band %d", b)
            lines!(ax, r.x, r.E_sunny[b, :]; linewidth=2, label=lab)
        end
        axislegend(ax, position=:rb, framevisible=false)
    end

    fig_path = joinpath(out_fig_dir, "analytical_vs_sunny_clean_dispersion.png")
    save(fig_path, fig)

    println("Saved outputs")
    println("-------------")
    println("figure       = ", fig_path)
    println("long CSV     = ", long_csv)
    println("summary CSV  = ", summary_csv)
    println("geometry CSV = ", geom_csv)
    println()
    println("Interpretation tip:")
    println("  This script uses a P1 lattice and explicitly lists all three positive")
    println("  J1 and J2 bond representatives matching the analytical form factors.")

    return (; figure=fig_path, long_csv, summary_csv, geometry_csv=geom_csv, results)
end

main()
