# Figure 5 ALAN Reallocation Reproducible Outputs

This folder contains a self-contained reproduction package for Figure 5:

- `policy_curve_gain_visible`: before-and-after female PA intensity and female-enabled footprint along the ALAN radiance gradient.
- `selected_city_action_change_map`: policy action map for the nine selected cities.

## Contents

- `data/policy_curve_profiles.csv`: binned curve data used for the ALAN response visualization.
- `data/city_action_units.csv`: selected-city spatial-unit data used for the action map.
- `data/city_action_summary.csv`: city-by-action summary table for checking map composition.
- `scripts/fig5_ALAN_reallocation.py`: plotting script.
- `figures/`: PNG and SVG outputs.

## Run

From this folder, run:

```bash
python scripts/fig5_ALAN_reallocation.py
```

The script reads the CSV files in `data/` and writes PNG and SVG figures to `figures/`.

## Outputs

- `figures/policy_curve_gain_visible.png`
- `figures/policy_curve_gain_visible.svg`
- `figures/selected_city_action_change_map.png`
- `figures/selected_city_action_change_map.svg`

The selected cities are Hong Kong, Dubai, Sydney, Auckland, Kuala Lumpur, Doha, Bangkok, Tehran and Dhaka.
