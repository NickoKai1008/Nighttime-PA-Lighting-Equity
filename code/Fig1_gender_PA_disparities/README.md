# Gender footprint figures

This folder contains a clean, reproducible package for selected gender-footprint figures.

## Files

- `scripts/gender_footprint.py`: standalone plotting script.
- `data/gender_grid_month_categories_200m.csv`: 200-m grid category input for Hong Kong and Dhaka.
- `data/daily_nighttime_footprint_area_by_city_sex.csv`: daily nighttime footprint area input for Hong Kong and Dhaka.
- `data/city_day_night_female_footprint_share_25cities.csv`: city-level summary for the 25-city footprint-share scatter.
- `maps/`: regenerated Hong Kong and Dhaka grid maps.
- `figures/`: regenerated density and 25-city footprint-share scatter figures.

## Run

```powershell
.\Fig1_gender_PA_disparities\scripts\gender_footprint.py
```

The script reads only local files from `data/` and writes to `figures/` and `maps/`.
