using TOML

"""
    load_toml_config(path)

Load a TOML configuration file and return the parsed dictionary.
"""
function load_toml_config(path::AbstractString)
    if !isfile(path)
        error("Could not find TOML config file: $path")
    end

    return TOML.parsefile(path)
end

"""
    load_best_fit_parameters(path)

Load the best-fit parameter TOML file used for YbZn2GaO5 neutron/magnetization modeling.
"""
function load_best_fit_parameters(path::AbstractString)
    return load_toml_config(path)
end

"""
    load_cofit_controls(path)

Load the co-fit run controls, such as data selection, weights, sampling,
plotting, and output options.
"""
function load_cofit_controls(path::AbstractString)
    return load_toml_config(path)
end

"""
    toml_symbol(value)

Convert string-valued TOML options such as `"tail_bgsub"` into Julia Symbols.
"""
toml_symbol(value::AbstractString) = Symbol(value)
toml_symbol(value::Symbol) = value