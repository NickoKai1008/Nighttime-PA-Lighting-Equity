# Figure 3 — Nighttime activity and road radiance

From this folder, run

```sh
python -m pip install -r requirements.txt
python plot.py
```

The script writes `output/fig3.png` (450 dpi), `output/fig3.svg` and `output/fig3.pdf`.

| Panel | Input in `data/` | Contents |
| --- | --- | --- |
| a, b | `pooled_curves.csv`, `bootstrap_bands.csv`, `female_thresholds.csv`, `grid_day_histograms.csv` | Pooled and land-use-specific fitted responses, pointwise 95% spatial-block bootstrap intervals, women’s flattening onset and L90, and grid-day counts per radiance bin. |
| c | `onset_ratios.csv` | Median female-to-male onset ratios and 95% intervals from paired spatial-block bootstrap draws. |
| d | `exposure_shares.csv` | Fractions below onset, between onset and L90, and above L90, weighted by grid-day count or by summed radiance. |
| e | `city_curves.csv`, `city_thresholds.csv` | Fitted responses, model-based pointwise 95% intervals and women’s local thresholds for the nine cities shown. |

The curve tables contain log relative activity (`f`) and log-scale interval bounds (`lo`, `hi`). The script exponentiates these values for display. Reference radiance is the group median in a–b and the city median in e. Radiance values in the CSVs are in W m⁻² sr⁻¹ μm⁻¹. The horizontal axes display radiance multiplied by 1,000.

Panel a computes the inset bars from adjacent differences in the fitted activity ratios at radiance values 2, 4, …, 16 in display units. The three photographs are in `images/`. Panel e uses the nighttime shared smooth (basis dimension 6) with city-specific deviations (basis dimension 5). The horizontal range is 0.08–12 in display units. Arrows mark curves or thresholds extending beyond the displayed range.
