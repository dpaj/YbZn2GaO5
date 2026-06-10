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


# ---------------------------------------------------------------------------
# Best-fit parameter helpers
# ---------------------------------------------------------------------------

function _lookup_toml_value(config::Dict, sections::Vector{String}, key::String)
    for section in sections
        if haskey(config, section)
            table = config[section]
            if table isa Dict && haskey(table, key)
                return table[key]
            end
        end
    end

    tried = join(["[$s].$key" for s in sections], ", ")
    error("Missing required TOML value. Tried: $tried")
end

function _lookup_positive_value_or_log10(
    config::Dict;
    sections::Vector{String},
    value_key::String,
    log10_key::String,
)
    for section in sections
        if !haskey(config, section)
            continue
        end

        table = config[section]
        if !(table isa Dict)
            continue
        end

        if haskey(table, value_key)
            return Float64(table[value_key])
        end

        if haskey(table, log10_key)
            return 10.0 ^ Float64(table[log10_key])
        end
    end

    tried_value = join(["[$s].$value_key" for s in sections], ", ")
    tried_log10 = join(["[$s].$log10_key" for s in sections], ", ")
    error("Missing required TOML value. Tried: $tried_value or $tried_log10")
end

"""
    cofit_initial_guess_kwargs(config)

Convert `configs/best_fit_parameters.toml` into keyword arguments for
`cofit_default_param_specs`.

This intentionally returns keyword names such as `initial_gzz`,
`initial_J1_meV`, etc., because those are the names expected by the
legacy co-fit parameter-spec constructor.
"""
function cofit_initial_guess_kwargs(config::Dict)
    physical_sections = ["physical", "initial_guess", "parameters"]

    extrinsic_sections = [
        "neutron_extrinsic",
        "extrinsic",
        "physical",
        "initial_guess",
        "parameters",
    ]

    return (;
        initial_gzz = Float64(_lookup_toml_value(config, physical_sections, "gzz")),
        initial_J1_meV = Float64(_lookup_toml_value(config, physical_sections, "J1_meV")),
        initial_J2_meV = Float64(_lookup_toml_value(config, physical_sections, "J2_meV")),
        initial_sigma_gzz = Float64(_lookup_toml_value(config, physical_sections, "sigma_gzz")),
        initial_sigma_J = Float64(_lookup_toml_value(config, physical_sections, "sigma_J")),

        initial_gzz2 = Float64(_lookup_toml_value(config, physical_sections, "gzz2")),
        initial_sigma_gzz2 = Float64(_lookup_toml_value(config, physical_sections, "sigma_gzz2")),

        initial_gperp_ratio = Float64(_lookup_toml_value(config, physical_sections, "gperp_ratio")),
        initial_chi_vv_muB_per_T = Float64(_lookup_toml_value(config, physical_sections, "chi_vv_muB_per_T")),

        initial_shared_r2 = _lookup_positive_value_or_log10(
            config;
            sections = extrinsic_sections,
            value_key = "second_kernel_relative_intensity",
            log10_key = "log10_second_kernel_relative_intensity",
        ),

        initial_neutron_global_scale = _lookup_positive_value_or_log10(
            config;
            sections = extrinsic_sections,
            value_key = "neutron_global_scale",
            log10_key = "log10_neutron_scale",
        ),
    )
end

"""
    print_initial_guess_kwargs(kwargs)

Print the initial guesses that will be passed to the co-fit parameter specs.
"""
function print_initial_guess_kwargs(kwargs::NamedTuple)
    println("Initial guesses loaded from TOML")
    println("--------------------------------")
    for name in propertynames(kwargs)
        println(rpad(String(name), 44), " = ", getproperty(kwargs, name))
    end
    println()
end