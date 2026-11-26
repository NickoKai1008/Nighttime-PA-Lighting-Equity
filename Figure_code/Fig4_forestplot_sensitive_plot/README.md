# Figure 4

Adjusted associations with nighttime physical activity and sensitivity analyses across 25 cities. Panel letters follow the manuscript.

## Draw the figure

With Python, Matplotlib and pandas installed, run from this directory.

```sh
python plot.py
```

The script writes `output/fig4.png` (300 dpi), `output/fig4.svg` and `output/fig4.pdf`.

| Panel | Input in `data/` | Quantity |
| --- | --- | --- |
| a | `coefficients.csv` | Sex-specific adjusted log coefficients, standard errors, P values and model-based 95% intervals |
| d | `night_day_gain.csv` | Night-to-day ratio of the fitted P99-to-P1 activity gain, with paired spatial-block bootstrap intervals |
| e | `clock_solar.csv` | Night-to-day gain ratios under clock-time and solar-elevation definitions |
| f | `perception_adjustment.csv` | Relative curve distance after adjusting for safety alone or safety and liveliness, shown as a percentage |
| g | `safety_terciles.csv` | Nighttime gain in the lowest versus highest within-city safety tercile, with paired spatial-block bootstrap intervals |

## Refit panel a

The refit was tested with R 4.5.3 and `mgcv` 1.9-4. Each sex-specific model uses 2,800,643 grid-days. Run from this directory.

```sh
Rscript model/fit.R
python plot.py --refit
```

The R script fits a quasi-Poisson model with a log link separately for women and men. It uses a linear term for log road radiance, standardized grid, weather and city covariates, a weekend indicator, a cubic regression spline for date and a city random intercept. It compares the 62 refitted coefficients and their uncertainty with `data/coefficients.csv` and writes `output/refit_coefficients.csv`. The second command uses these estimates for panel a.

The model inputs are in `model/data/`. Each row of `observations.csv.gz` gives female and male nighttime activity for one grid-day. Within each city, `cell` links the observation to `cells.csv.gz`, and `period` links it to `periods.csv.gz`. `cities.csv` supplies city-level covariates. `scaling.csv` supplies the centers and standard deviations for grid and weather covariates. City covariates are standardized across the 25 cities.
