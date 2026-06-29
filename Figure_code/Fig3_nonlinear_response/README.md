# Fig. 3 nonlinear ALAN response plotting package

This folder contains figure-ready response-curve tables and a standalone plotting script for redrawing the Fig. 3 nonlinear ALAN response panels.

## Contents

- `data/panel_a_pooled_response.csv`: pooled female nighttime PA response curve.
- `data/panel_a_pooled_thresholds.csv`: pooled threshold table used for vertical guide lines.
- `data/panel_b_city_conditioned_response_9cities.csv`: city-conditioned response curves for the nine displayed cities.
- `data/panel_b_city_thresholds_9cities.csv`: city-conditioned threshold table used for vertical guide lines.
- `data/panel_c_landuse_response.csv`: functional-context response curves.
- `data/panel_c_landuse_thresholds.csv`: functional-context threshold table used for vertical guide lines.
- `data/observed_ntl_histograms_female_dropzero_log_raw.csv`: observed NTL histograms used in panels a and c.
- `scripts/fig3_nonlinear_response.py`: standalone plotting script.


## Run

```powershell
.\scripts\fig3_nonlinear_response.py
```

The script writes PNG and SVG outputs to `figures/`.
