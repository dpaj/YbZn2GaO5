module YZGOCofit

include("parameters.jl")

# Feature extraction is optional for lightweight validation/replay scripts.
# If the file is present in the repo, include and export its public helpers.
if isfile(joinpath(@__DIR__, "feature_extraction.jl"))
    include("feature_extraction.jl")
end

export load_toml_config
export load_best_fit_parameters
export load_cofit_controls
export toml_symbol
export cofit_initial_guess_kwargs
export print_initial_guess_kwargs
export load_plot_2d_controls
export canonical_model_parameters
export load_canonical_model_parameters
export canonical_model_parameters_dict
export print_canonical_model_parameters

# Optional feature-extraction exports. These symbols are defined when
# `src/feature_extraction.jl` is present.
export run_neutron_magnetization_feature_extraction
export fe_load_neutron_scans
export fe_fit_neutron_features_for_field
export fe_fit_magnetization_features

end
