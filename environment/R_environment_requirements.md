# R analysis requirements

The restricted-data analysis scripts in `code/00_data_preparation` through
`code/07_supplementary` use R and the following packages:

- `data.table`, `tidyverse`, `dplyr`, `tidyr`, `readr`
- `mgcv` 1.9-1 and `gratia`
- `lme4`, `lmerTest`, `broom.mixed`, `performance`
- `ggplot2`, `patchwork`, `scales`, `sf`
- `openxlsx`, `writexl`, `yaml`

The `mgcv` version is the version reported for the fitted GAMMs in the
manuscript. The remaining packages are listed because they are imported by the
released scripts; exact versions were not archived. The public figure demo is
separate and uses `environment/environment.yml`.
