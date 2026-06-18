#!/usr/bin/env julia

# Standalone Sunny.jl KPM benchmark for the YbZn2GaO5 field-polarized
# dispersive component.
#
# Purpose: isolate Sunny's core KPM workflow without repository TOMLs,
# experimental data loading, background subtraction, scaling, or plotting.

using Sunny
using LinearAlgebra
using Random
using Printf
using Statistics
using Dates

# -----------------------------------------------------------------------------
# User knobs
# -----------------------------------------------------------------------------

const SEED = 20260618
const OUTPUT_DIR = "sunny_kpm_core_benchmark_output"

const A_ANGSTROM = 3.376
const C_ANGSTROM = 21.96
const SPIN_S = 0.5

# Canonical YZGO best-fit values. Multipliers let this script be used for
# stress tests without reading project TOMLs.
const GZZ = 3.784584444
const J1_MEV = 0.2299824987
const J2_MEV = 0.008722296371
const SIGMA_GZZ = 0.3859237343
const SIGMA_J = 0.2396873058
const SIGMA_GZZ_MULTIPLIER = 1.0
const SIGMA_J_MULTIPLIER = 1.0

# Nonzero transverse g is needed for the neutron transverse response in Sunny's
# ssf_perp measure. This is an intensity gauge for benchmarking, not a fitted
# physical g_perp.
const SUNNY_TRANSVERSE_GXY = 1.0

const FIELD_T = 9.0
const FIELD_DIRECTION = [0.0, 0.0, 1.0]

# Build a 3x3 magnetic cell and repeat it to a 36x36 inhomogeneous system.
const DIMS = (3, 3, 1)
const REPEAT_FACTOR = (12, 12, 1)

const INITIALIZE_FIELD_POLARIZED = true
const RELAX_GROUND_STATE = true
const MAXITERS = 5000

const KPM_TOL = 0.05
const KERNEL_FWHM_MEV = 0.08
const ENERGY_MIN_MEV = 0.6
const ENERGY_MAX_MEV = 2.8
const N_ENERGY = 161

# Q grid: measured cut-volume grid crossed with a small Gaussian-like momentum
# resolution grid. L is intentionally one-point here; add L later for intensity
# integration once the core timing path is validated.
const Q_CENTER = [0.0, 1.0, 0.0]
const MEASURED_HALF_WIDTH = [0.05, 0.05, 0.0]
const RESOLUTION_SIGMA = [0.035, 0.035, 0.0]
const RESOLUTION_NSIGMA = 2.0
const N_MEASURED = (5, 5, 1)
const N_RESOLUTION = (5, 5, 1)

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

struct TimerRows
    rows::Vector{NamedTuple}
end

TimerRows() = TimerRows(NamedTuple[])

function timed!(f, timers::TimerRows, stage::AbstractString; q_samples::Int=0, note::AbstractString="")
    t0 = time_ns()
    result = f()
    elapsed = (time_ns() - t0) / 1e9
    push!(timers.rows, (; stage=String(stage), seconds=elapsed, q_samples=q_samples, note=String(note)))
    @printf("%-24s %10.3f s", stage, elapsed)
    q_samples > 0 && @printf("   (%d Q, %.6f s/Q)", q_samples, elapsed / q_samples)
    isempty(note) || @printf("   %s", note)
    println()
    return result
end

function write_profile_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "stage,seconds,q_samples,seconds_per_q,note")
        for r in rows
            spq = r.q_samples > 0 ? r.seconds / r.q_samples : NaN
            note = replace(r.note, '"' => "'", ',' => ';')
            @printf(io, "%s,%.9g,%d,%.9g,%s\n", r.stage, r.seconds, r.q_samples, spq, note)
        end
    end
    return path
end

function write_spectrum_csv(path::AbstractString, energies, spectrum)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "energy_meV,intensity_avg")
        for i in eachindex(energies)
            @printf(io, "%.10g,%.10g\n", Float64(energies[i]), Float64(spectrum[i]))
        end
    end
    return path
end

function write_qgrid_csv(path::AbstractString, qpts, weights)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "H,K,L,weight")
        for i in eachindex(qpts)
            q = qpts[i]
            @printf(io, "%.10g,%.10g,%.10g,%.10g\n", q[1], q[2], q[3], weights[i])
        end
    end
    return path
end

function write_summary_txt(path::AbstractString, profile_rows, spectrum; extra="")
    mkpath(dirname(path))
    total = sum(r.seconds for r in profile_rows)
    intensity_rows = filter(r -> r.stage == "kpm_intensities", profile_rows)
    kpm_time = isempty(intensity_rows) ? NaN : sum(r.seconds for r in intensity_rows)
    open(path, "w") do io
        println(io, "Sunny KPM core benchmark")
        println(io, "timestamp = ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
        println(io, "Sunny version = ", try string(pkgversion(Sunny)) catch; "unknown" end)
        println(io)
        println(io, "field_T = ", FIELD_T)
        println(io, "dims = ", DIMS)
        println(io, "repeat_factor = ", REPEAT_FACTOR)
        println(io, "system_size = ", (DIMS[1]*REPEAT_FACTOR[1], DIMS[2]*REPEAT_FACTOR[2], DIMS[3]*REPEAT_FACTOR[3]))
        println(io, "sigma_J_used = ", SIGMA_J_MULTIPLIER * SIGMA_J)
        println(io, "sigma_gzz_used = ", SIGMA_GZZ_MULTIPLIER * SIGMA_GZZ)
        println(io, "q_points = ", N_MEASURED[1]*N_MEASURED[2]*N_MEASURED[3]*N_RESOLUTION[1]*N_RESOLUTION[2]*N_RESOLUTION[3])
        println(io, "n_energy = ", N_ENERGY)
        println(io, "kpm_tol = ", KPM_TOL)
        println(io, "kernel_fwhm_meV = ", KERNEL_FWHM_MEV)
        println(io)
        println(io, "total_timed_seconds = ", total)
        println(io, "kpm_intensities_seconds = ", kpm_time)
        println(io, "kpm_fraction = ", isfinite(kpm_time) ? kpm_time / total : NaN)
        println(io, "spectrum_sum = ", sum(spectrum))
        isempty(extra) || println(io, extra)
    end
    return path
end

# Positive representative shell offsets for the triangular P1 cell, matched to
# the analytical YZGO convention:
# Δ1 = 6 - 2[cos(2πH) + cos(2πK) + cos(2π(H+K))]
# Δ2 = 6 - 2[cos(2π(H-K)) + cos(2π(2H+K)) + cos(2π(H+2K))]
j1_shell_offsets() = ([1, 0, 0], [0, 1, 0], [1, 1, 0])
j2_shell_offsets() = ([1, -1, 0], [2, 1, 0], [1, 2, 0])

function field_unit_vector()
    u = Float64.(FIELD_DIRECTION)
    n = norm(u)
    n > 0 || error("FIELD_DIRECTION must be nonzero")
    return u ./ n
end

function set_g_tensor_all_sites!(sys; gzz::Real, gxy::Real)
    for site in eachsite(sys)
        sys.gs[site] = [Float64(gxy) 0.0 0.0;
                        0.0 Float64(gxy) 0.0;
                        0.0 0.0 Float64(gzz)]
    end
    return sys
end

function build_clean_system()
    units = Units(:meV, :angstrom)
    latvecs = lattice_vectors(A_ANGSTROM, A_ANGSTROM, C_ANGSTROM, 90, 90, 120)
    cryst = Crystal(latvecs, [[0.0, 0.0, 0.0]], 1; types=["Yb"])
    sys = System(cryst, [1 => Moment(s=SPIN_S, g=1.0)], :dipole; dims=DIMS, seed=SEED)

    set_g_tensor_all_sites!(sys; gzz=GZZ, gxy=SUNNY_TRANSVERSE_GXY)
    for off in j1_shell_offsets()
        set_exchange!(sys, J1_MEV, Bond(1, 1, off))
    end
    for off in j2_shell_offsets()
        set_exchange!(sys, J2_MEV, Bond(1, 1, off))
    end

    set_field!(sys, collect(field_unit_vector() .* (FIELD_T * units.T)))
    return (; sys, cryst, units)
end

function make_inhomogeneous_with_disorder!(sys)
    sys = to_inhomogeneous(repeat_periodically(sys, REPEAT_FACTOR))
    rng = MersenneTwister(SEED + 7919)

    sigma_J_used = SIGMA_J_MULTIPLIER * SIGMA_J
    for off in j1_shell_offsets()
        for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
            Jij = J1_MEV * (1.0 + sigma_J_used * randn(rng))
            set_exchange_at!(sys, Jij, s1, s2; offset=o)
        end
    end
    for off in j2_shell_offsets()
        for (s1, s2, o) in symmetry_equivalent_bonds(sys, Bond(1, 1, off))
            Jij = J2_MEV * (1.0 + sigma_J_used * randn(rng))
            set_exchange_at!(sys, Jij, s1, s2; offset=o)
        end
    end

    sigma_gzz_used = SIGMA_GZZ_MULTIPLIER * SIGMA_GZZ
    for site in eachsite(sys)
        gi = max(0.0, GZZ + sigma_gzz_used * randn(rng))
        sys.gs[site] = [SUNNY_TRANSVERSE_GXY 0.0 0.0;
                        0.0 SUNNY_TRANSVERSE_GXY 0.0;
                        0.0 0.0 gi]
    end
    return sys
end

function initialize_field_polarized!(sys)
    u = field_unit_vector()
    for site in eachsite(sys)
        set_dipole!(sys, u, site)
    end
    return sys
end

function linspace_offsets(half_width::Real, n::Int)
    n <= 1 && return [0.0]
    return collect(range(-Float64(half_width), Float64(half_width); length=n))
end

function gaussian_offsets(sigma::Real, nsigma::Real, n::Int)
    if sigma <= 0 || n <= 1 || !isfinite(sigma)
        return [0.0], [1.0]
    end
    xs = collect(range(-Float64(nsigma)*Float64(sigma), Float64(nsigma)*Float64(sigma); length=n))
    ws = exp.(-0.5 .* (xs ./ Float64(sigma)).^2)
    ws ./= sum(ws)
    return xs, ws
end

function build_q_grid()
    mh = linspace_offsets(MEASURED_HALF_WIDTH[1], N_MEASURED[1])
    mk = linspace_offsets(MEASURED_HALF_WIDTH[2], N_MEASURED[2])
    ml = linspace_offsets(MEASURED_HALF_WIDTH[3], N_MEASURED[3])

    rh, wh = gaussian_offsets(RESOLUTION_SIGMA[1], RESOLUTION_NSIGMA, N_RESOLUTION[1])
    rk, wk = gaussian_offsets(RESOLUTION_SIGMA[2], RESOLUTION_NSIGMA, N_RESOLUTION[2])
    rl, wl = gaussian_offsets(RESOLUTION_SIGMA[3], RESOLUTION_NSIGMA, N_RESOLUTION[3])

    qpts = Vector{Vector{Float64}}()
    weights = Float64[]
    measured_weight = 1.0 / (length(mh) * length(mk) * length(ml))

    for dhm in mh, dkm in mk, dlm in ml
        for (i, dhr) in pairs(rh), (j, dkr) in pairs(rk), (k, dlr) in pairs(rl)
            push!(qpts, [Q_CENTER[1] + dhm + dhr,
                         Q_CENTER[2] + dkm + dkr,
                         Q_CENTER[3] + dlm + dlr])
            push!(weights, measured_weight * wh[i] * wk[j] * wl[k])
        end
    end
    weights ./= sum(weights)
    return qpts, weights
end

function intensity_data_matrix(Iraw)
    if Iraw isa AbstractMatrix
        return Matrix(Iraw)
    elseif hasproperty(Iraw, :data)
        return Matrix(getproperty(Iraw, :data))
    elseif hasproperty(Iraw, :intensities)
        return Matrix(getproperty(Iraw, :intensities))
    else
        error("Could not convert intensities output of type $(typeof(Iraw)) to a matrix")
    end
end

function normalize_intensity_matrix(Iraw, nE::Int, nQ::Int)
    I = intensity_data_matrix(Iraw)
    if size(I) == (nE, nQ)
        return I
    elseif size(I) == (nQ, nE)
        return transpose(I) |> Matrix
    else
        error("Unexpected intensities matrix size $(size(I)); expected ($nE,$nQ) or ($nQ,$nE)")
    end
end

# -----------------------------------------------------------------------------
# Main benchmark
# -----------------------------------------------------------------------------

function main()
    timers = TimerRows()
    mkpath(OUTPUT_DIR)

    @printf("Sunny KPM core benchmark\n")
    @printf("------------------------\n")
    @printf("Sunny version       = %s\n", try string(pkgversion(Sunny)) catch; "unknown" end)
    @printf("field_T             = %.3f\n", FIELD_T)
    @printf("system_size         = (%d, %d, %d)\n", DIMS[1]*REPEAT_FACTOR[1], DIMS[2]*REPEAT_FACTOR[2], DIMS[3]*REPEAT_FACTOR[3])
    @printf("sigma_J used        = %.9g\n", SIGMA_J_MULTIPLIER * SIGMA_J)
    @printf("sigma_gzz used      = %.9g\n", SIGMA_GZZ_MULTIPLIER * SIGMA_GZZ)
    @printf("grid                = measured %dx%dx%d, resolution %dx%dx%d\n",
            N_MEASURED..., N_RESOLUTION...)
    println()

    built = timed!(timers, "build_clean_system") do
        build_clean_system()
    end

    sys = timed!(timers, "enlarge_and_disorder") do
        make_inhomogeneous_with_disorder!(built.sys)
    end

    timed!(timers, "initialize_spins") do
        INITIALIZE_FIELD_POLARIZED ? initialize_field_polarized!(sys) : randomize_spins!(sys)
    end

    if RELAX_GROUND_STATE
        timed!(timers, "minimize_energy") do
            minimize_energy!(sys; maxiters=MAXITERS)
        end
    end

    swt = timed!(timers, "construct_kpm") do
        SpinWaveTheoryKPM(sys; measure=ssf_perp(sys), tol=KPM_TOL)
    end

    energies = collect(range(ENERGY_MIN_MEV, ENERGY_MAX_MEV; length=N_ENERGY))
    kernel = gaussian(fwhm=KERNEL_FWHM_MEV)

    qpts, weights = timed!(timers, "build_q_grid") do
        build_q_grid()
    end
    nQ = length(qpts)

    Iraw = timed!(timers, "kpm_intensities"; q_samples=nQ, note="CPU Sunny") do
        intensities(swt, qpts; energies=energies, kernel=kernel, kT=0.0, verbose=false)
    end

    spectrum = timed!(timers, "average_spectrum"; q_samples=nQ) do
        I = normalize_intensity_matrix(Iraw, length(energies), nQ)
        I * weights
    end

    profile_path = joinpath(OUTPUT_DIR, "profile.csv")
    spectrum_path = joinpath(OUTPUT_DIR, "spectrum.csv")
    qgrid_path = joinpath(OUTPUT_DIR, "qgrid.csv")
    summary_path = joinpath(OUTPUT_DIR, "summary.txt")

    timed!(timers, "write_outputs") do
        write_profile_csv(profile_path, timers.rows)
        write_spectrum_csv(spectrum_path, energies, spectrum)
        write_qgrid_csv(qgrid_path, qpts, weights)
        write_summary_txt(summary_path, timers.rows, spectrum)
    end

    println()
    println("Saved:")
    println("  ", profile_path)
    println("  ", spectrum_path)
    println("  ", qgrid_path)
    println("  ", summary_path)
    return nothing
end

main()
