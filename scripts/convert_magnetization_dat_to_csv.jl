#!/usr/bin/env julia
# Convert Quantum Design .dat magnetometry into the two-column CSV the repo's loader expects.
# DATA ONLY, no model, seconds to run.
#
#   julia --project=. scripts/convert_magnetization_dat_to_csv.jl
#
# WHY A SCRIPT RATHER THAN A LOADER
#
# Teaching `sv_read_magnetization_csv` to parse QD files would bury the two choices that matter
# inside library code. Both are choices, not conventions, and both have already caused a wrong
# answer, so they belong in a visible, re-runnable, version-controlled conversion step:
#
#   1. WHICH MOMENT COLUMN. In the MPMS3 file `Moment (emu)` is EMPTY. The moment lives in
#      `DC Moment Fixed Ctr` (sample held at the nominal position) and `DC Moment Free Ctr`
#      (position floated as a fit parameter). The sample sat 3.06 mm off centre, so the fixed fit
#      solved at the wrong place and under-read by a constant 1.486. The file's own fit-quality
#      columns settle it: fixed 0.268, free 0.957. FREE is correct. The digitised curve used until
#      now came from FIXED, which is the entire origin of the long-standing 1.5x scale problem.
#
#   2. WHICH TEMPERATURE COLUMN. In the same file `Temperature (K)` reads 1.56 K -- the MPMS3
#      CHAMBER. Sample temperature is in a column named `He3 temp`, reading 0.42-0.50 K. Using the
#      obvious column would be wrong by a factor of four.
#
# The DynaCool VSM files have neither problem: `Moment (emu)` is populated and
# `Temperature (K)` is the sample.
#
# UNITS. QD reports emu. mu_B/Yb = emu / (5585 * moles), moles = mass_g / 453.53, using
# M(YbZn2GaO5) = 453.53 g/mol (Yb 173.05 + 2xZn 130.76 + Ga 69.72 + 5xO 80.00), one Yb per formula
# unit, and 5585 emu/mol = N_A * mu_B. Mass comes from the FILENAME -- the SAMPLE_MASS header field
# is blank in every one of these files.
#
# The output keeps B_T and M_muB_per_Yb as columns 1 and 2 so `sv_read_magnetization_csv` reads it
# unchanged, and carries the provenance alongside. For the MPMS3 file the superseded fixed-centring
# value is preserved as an extra column, so the correction stays auditable from the data itself.

using Printf, Statistics

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const MOL = 453.53      # g/mol, YbZn2GaO5
const NA_MUB = 5585.0   # emu/mol per mu_B per formula unit

"Read a QD .dat file: Dict of column name => Vector{Float64}, NaN where non-numeric."
function read_qd(path)
    lines = readlines(path)
    i = findfirst(l -> strip(l) == "[Data]", lines)
    i === nothing && error("No [Data] section in $path")
    hdr = strip.(String.(split(lines[i+1], ',')))
    cols = Dict{String,Vector{Float64}}(h => Float64[] for h in hdr if !isempty(h))
    for l in lines[i+2:end]
        f = String.(split(l, ','))
        length(f) < length(hdr) && continue
        for (k, h) in enumerate(hdr)
            isempty(h) && continue
            v = tryparse(Float64, strip(f[k]))
            push!(cols[h], v === nothing ? NaN : v)
        end
    end
    return cols
end

emu_to_muB(mass_mg) = 1.0 / (NA_MUB * mass_mg * 1e-3 / MOL)

# file, mass (mg), moment column, error column, sample-temperature column, label, output name
const SPECS = [
    (rel = "data/magnetization/mpms3_0p4K/YbZnGaO_Bpara001_12.2mg_MvH_0.4K_10042025.dat",
     mass = 12.2, mom = "DC Moment Free Ctr (emu)", err = "DC Moment Err Free Ctr (emu)",
     temp = "He3 temp", also = "DC Moment Fixed Ctr (emu)",
     label = "0.42 K, B || c, MPMS3 + He-3 insert, free centring",
     out = "YZGO_MvH_0p42K_Bparc_mpms3.csv"),
    (rel = "data/magnetization/ppms_2p5K/YZGO_BparaC_4.81MG_2.5K_06242026.DAT",
     mass = 4.81, mom = "Moment (emu)", err = "M. Std. Err. (emu)",
     temp = "Temperature (K)", also = nothing,
     label = "2.5 K, B || c, DynaCool VSM",
     out = "YZGO_MvH_2p5K_Bparc_vsm.csv"),
    (rel = "data/magnetization/ppms_2p5K/YZGO_BperpC_7.7MG_2.5K_06242026.DAT",
     mass = 7.7, mom = "Moment (emu)", err = "M. Std. Err. (emu)",
     temp = "Temperature (K)", also = nothing,
     label = "2.5 K, B perpendicular to c, DynaCool VSM",
     out = "YZGO_MvH_2p5K_Bperpc_vsm.csv"),
]

for sp in SPECS
    path = joinpath(REPO_ROOT, sp.rel)
    isfile(path) || (@warn "missing, skipping" file=sp.rel; continue)
    c = read_qd(path)
    haskey(c, sp.mom) || error("No column '$(sp.mom)' in $(sp.rel)")
    conv = emu_to_muB(sp.mass)
    B = c["Magnetic Field (Oe)"] ./ 1e4
    M = c[sp.mom] .* conv
    E = haskey(c, sp.err) ? c[sp.err] .* conv : fill(NaN, length(B))
    T = haskey(c, sp.temp) ? c[sp.temp] : fill(NaN, length(B))
    A = sp.also !== nothing && haskey(c, sp.also) ? c[sp.also] .* conv : nothing

    keep = isfinite.(B) .& isfinite.(M)
    outp = joinpath(REPO_ROOT, "data/magnetization", sp.out)
    open(outp, "w") do io
        # Columns 1 and 2 are what sv_read_magnetization_csv reads; the rest is provenance.
        extra = A === nothing ? "" : ",M_muB_per_Yb_fixed_ctr_SUPERSEDED"
        println(io, "B_T,M_muB_per_Yb,M_err_muB_per_Yb,T_sample_K,moment_column,mass_mg,source_file" * extra)
        for k in eachindex(B)
            keep[k] || continue
            @printf(io, "%.6g,%.6g,%.6g,%.4f,%s,%.2f,%s", B[k], M[k], E[k], T[k],
                    replace(sp.mom, "," => ";"), sp.mass, basename(sp.rel))
            A === nothing || @printf(io, ",%.6g", A[k])
            println(io)
        end
    end

    at(t) = M[argmin(abs.(B .- t))]
    @printf("%-46s -> %s\n", sp.label, sp.out)
    @printf("   %d pts, T = %.3f-%.3f K, B = %.2f-%.2f T, M = %.4f-%.4f uB/Yb",
            count(keep), minimum(filter(isfinite, T)), maximum(filter(isfinite, T)),
            minimum(B), maximum(B), minimum(M[keep]), maximum(M[keep]))
    @printf("\n   M(6.975 T) = %.4f uB/Yb\n", at(6.975))
    if A !== nothing
        r = filter(isfinite, (M ./ A)[A .> 0.02])
        @printf("   superseded fixed-centring column is a factor %.4f +- %.4f lower\n",
                mean(r), std(r))
    end
end

println()
println("NOTE the configs still point at data/magnetization/YZGO_MvB_black_curve_digitized_visible.csv,")
println("which is superseded -- digitised from a figure AND from the wrong moment column. Switching")
println("them is deliberate and should follow a digitised-versus-primary comparison, because it")
println("changes every absolute M(H) number in the repo (shape fits are unaffected: A_M absorbs it).")
