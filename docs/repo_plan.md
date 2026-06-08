# Repo plan

Near-term goal: make a clean, shareable code base for the YbZn2GaO5 neutron + magnetization co-fit.

## Current priorities

1. Preserve the existing working co-fit analysis.
2. Move initial guesses and extrinsic controls into visible config files.
3. Separate reusable model/fitting functions from runnable scripts.
4. Add 2D data-versus-model plotting.
5. Add comparison of 3.32 meV and 4.65 meV data/background handling.
6. Add feature extraction as an interpretable intermediate analysis layer.
7. Add Sunny.jl final-parameter comparison as an independent validation calculation.

## Modeling strategy

The optimizer uses a fast analytical field-polarized model. Sunny.jl will be included as a slower validation backend for the final fitted model parameters.