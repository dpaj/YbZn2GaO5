using TOML

"""
    load_best_fit_parameters(path)

Load the best-fit parameter TOML file used for YbZn2GaO5 neutron/magnetization modeling.

This keeps the scientific/fitting parameters outside the scripts so that:
  - initial guesses are easy to find,
  - plotting scripts and fitting scripts use the same numbers,
  - future Sunny and synthetic-data scripts can use the same parameter source.
"""
function load_best_fit_parameters(path::AbstractString)
    if !isfile(path)
        error("Could not find parameter file: $path")
    end

    return TOML.parsefile(path)
end