# Data dictionary

        Key fields:

        - `female_sum_night`, `male_sum_night`: synthetic demo PA counts in the nighttime window; restricted in the real analysis.
        - `log_female_night`, `log_male_night`: log-transformed nighttime PA volume used for exposure-response modelling on positive grid-days.
        - `NTL_mean_grid`: synthetic demo ALAN radiance field; restricted in the real analysis.
        - `NDVI`, `built_density_2020`, `ssx_choice_r5000`, `walkability_*`: non-restricted built-environment covariates.
        - `land_use_type`: functional context used to estimate response-supported ALAN ranges.
        - `female_to_male_ratio_change_pct`: post-hoc city-level gender-balance diagnostic from summed female PA divided by summed male PA after applying simulated ALAN changes.
