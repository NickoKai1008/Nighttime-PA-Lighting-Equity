# Models for nighttime light and physical activity

The models use road-weighted nighttime light (NTL) radiance and daily counts of recorded physical activity (PA). Each fit uses one sex and one period. The examples below reproduce the model specifications for nighttime activity across all 25 cities.

| Analysis | R script and model option | Specification used in the paper |
| --- | --- | --- |
| Nonlinear NTL response | [fit_quasipoisson_shared_ntl_curve.R](fit_quasipoisson_shared_ntl_curve.R), `shared` | Main curves (Fig. 3a,b). Quasi-Poisson with a log link, a cubic regression spline of log radiance (`k = 6`) and city fixed intercepts. |
| Nonlinear NTL response by city | [fit_quasipoisson_city_deviation_curves.R](fit_quasipoisson_city_deviation_curves.R) | Fig. 3e and Supplementary Fig. S5b. Within each sex, a shared nighttime smooth (`k = 6`) and penalized city-specific deviations (`k_dev = 5`) are fitted jointly across all 25 cities. |
| Nonlinear NTL response with city covariates | [fit_multilevel_quasipoisson_linear_spline_and_two_part.R](fit_multilevel_quasipoisson_linear_spline_and_two_part.R), `FULL_EVEN` | Complementary multilevel curves (Supplementary Fig. S6a). A spline of log radiance (`k = 8`), city-level covariates and a city random intercept. |
| Linear NTL response | `fit_multilevel_quasipoisson_linear_spline_and_two_part.R`, `FULL_LIN` | Fig. 4a. The same multilevel quasi-Poisson model with a single coefficient for log radiance. |
| Two-part response | `fit_multilevel_quasipoisson_linear_spline_and_two_part.R`, `TP_OCC,TP_POS` | Main-text decomposition and Supplementary Method S12. A binomial logit model estimates the probability of a positive record. A Gaussian identity-link model estimates mean log(PA) among positive records. Both use the multilevel adjustment set and a log-radiance spline (`k = 8`). |

All specifications include a weekend indicator and a date smooth (`k = 5`). The count and occurrence models retain zero-record grid-days.

## Run the models

Use R 4.5.2 with `mgcv` 1.9-3, `data.table`, `jsonlite` and `openssl`. Set `V2_EX` to the prepared NDVI100 execution directory. The drivers look for these files relative to that directory:

```text
code/v2_model_common.R
code/v2_fit.R
cache/raw/all25__night__pooled__variant_expanded/RAW_WIDE.rds
cache/raw/all25__night__pooled__variant_expanded/RAW_WIDE.json
inputs/expanded/all25/night/METADATA.rds
inputs/FOLD_MAP.csv
```

The two R files in [code/](code/) provide the original model and data-preparation functions. For the multilevel models, `CITY` points to `17_city_level_12.csv`, the table of 12 city covariates keyed by `city`.

From this directory, run the following in a shell with `Rscript` on its path. Use a fresh output directory for each command.

```sh
export V2_EX="/path/to/NDVI100/replacement_20261109/execution"
CITY="/path/to/17_city_level_12.csv"

KNOTS=even Rscript --vanilla fit_quasipoisson_shared_ntl_curve.R all25 shared night results/shared_k6 6 female,male
for sex in female male; do
  KNOTS=even XSCALE=log BINCHECK=1 Rscript --vanilla fit_quasipoisson_city_deviation_curves.R all25 night "results/city_k6_$sex" 5 "$sex" 6
done
Rscript --vanilla fit_multilevel_quasipoisson_linear_spline_and_two_part.R all25 night 8 FULL_EVEN results/multilevel_spline "$CITY" female,male
Rscript --vanilla fit_multilevel_quasipoisson_linear_spline_and_two_part.R all25 night 8 FULL_LIN results/multilevel_linear "$CITY" female,male
Rscript --vanilla fit_multilevel_quasipoisson_linear_spline_and_two_part.R all25 night 8 TP_OCC,TP_POS results/two_part "$CITY" female,male
```

For the shared curve, `KNOTS=even` places `k − 1` knots evenly from the minimum to the 99th percentile of log radiance and the final knot at the maximum. `FULL_EVEN` applies this rule internally. The two-part fit uses `mgcv`’s default cubic-regression-spline knots. The `8` in the `FULL_LIN` command preserves the archived invocation; its NTL term has one coefficient.

In Supplementary Fig. S5a, each city panel shows the pooled curve over that city's observed radiance range. In Fig. 3e and Supplementary Fig. S5b, each curve combines the shared smooth with its city's deviation and is centred at the local median radiance. `CITY_CURVES_V2_*.csv` contains the city estimates and their model-based pointwise 95% intervals. Fig. 3e displays nine of the 25 fitted cities.

Each entry point writes a progress CSV with one row per fit. The shared and multilevel drivers write smooth-term estimates to `curves/`. The multilevel driver also writes `MULTILEVEL_CITY_COEF.csv`. For the complete linear forest plot, including the weekend coefficient, use the [Fig. 4 panel-a refit](../fig4/README.md#refit-panel-a).
