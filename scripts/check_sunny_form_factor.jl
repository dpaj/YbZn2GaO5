# scripts/check_sunny_form_factor.jl
#
# Quick diagnostic for the magnetic form-factor setting used by the Sunny KPM
# validation.  For source = sunny_builtin, it verifies Sunny FormFactor label
# construction and prints manual fallback values only as a diagnostic.
#
# Run from repo root:
#   julia --project=. scripts/check_sunny_form_factor.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "YZGOCofit.jl"))
using .YZGOCofit

include(joinpath(REPO_ROOT, "src", "sunny_validation.jl"))
using .SunnyValidation
using Printf

controls = SunnyValidation.sv_load_controls(REPO_ROOT)
ff = SunnyValidation.sv_form_factor_controls(controls)
lat = SunnyValidation.sv_form_factor_lattice_controls(controls)

println("Sunny neutron form-factor diagnostic")
println("-------------------------------------")
println("enabled       = ", ff.enabled)
println("source        = ", ff.source)
println("ion           = ", ff.ion)
println("candidate_ions= ", ff.candidate_ions)
println("manual_include_j2    = ", ff.include_j2)
println("manual_j2_coefficient= ", ff.c2)
println("manual_apply_as      = ", ff.apply_as)
println(@sprintf("manual lattice       = a %.6g A, c %.6g A, gamma %.6g deg", lat.a_A, lat.c_A, lat.gamma_deg))

if ff.enabled && ff.source == :sunny_builtin
    pairs = SunnyValidation.sv_builtin_formfactor_pairs(controls)
    println("Sunny built-in FormFactor pairs constructed successfully: ", pairs)
    println("The f(Q) columns below are the manual_yb3 fallback diagnostic only; the KPM intensity uses Sunny's built-in form factor internally.")
end
println()
println(rpad("point", 16), lpad("H", 10), lpad("K", 10), lpad("L", 10), lpad("|Q| A^-1", 14), lpad("manual f(Q)", 16), lpad("manual |f|^2", 16))

points = [
    ("Γ", [0.0, 0.0, 0.0]),
    ("K", [1/3, 1/3, 0.0]),
    ("M", [0.5, 0.0, 0.0]),
    ("Γ1", [1.0, 0.0, 0.0]),
    ("K1", [-1/3, 2/3, 0.0]),
]

for (name, q) in points
    Q = SunnyValidation.sv_qmag_Ainv(q, controls)
    f = SunnyValidation.sv_yb3_form_factor(Q; include_j2=ff.include_j2, c2=ff.c2)
    w = f^2
    println(rpad(name, 16), @sprintf("%10.5f%10.5f%10.5f%14.6f%16.6f%16.6f", q[1], q[2], q[3], Q, f, w))
end
