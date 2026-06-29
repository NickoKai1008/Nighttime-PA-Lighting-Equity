# S1 gender inequality quadrant figure

This folder is a minimal reproducible package for `Fig2_nighttime_gender_gap`.

## Files

- `scripts/plot_gender_inequality_quadrant.py`: standalone plotting script. It reads only the local CSV and outputs PNG and SVG.
- `figures/`: regenerated PNG and editable SVG.


- `data/lorenz_curve_points.csv`: city-level population-normalized Lorenz curve coordinates for nighttime female and male activity.
- `data/night_gini_by_city.csv`: city-level Gini values with GII group and geographic region labels.
- `scripts/Fig2c_lorenz_curves.py`: one plotting script that generates all figures.
- `figures/`: PNG and SVG outputs.


## Run

```powershell
.\Fig2_nighttime_gender_gap\scripts\plot_gender_inequality_quadrant.py
```


```bash
python scripts/Fig2c_lorenz_curves.py
```
