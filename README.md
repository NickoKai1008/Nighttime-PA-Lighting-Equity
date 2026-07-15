# Nighttime lighting and gender equity in physical activity

This repository provides public, release-safe materials for reproducing summary tables, diagnostics and demo outputs from a study of nighttime lighting and gendered physical activity.

Restricted activity records and satellite imagery are not redistributed. Public materials include non-restricted derived covariates, aggregate result tables, data dictionaries and synthetic demo inputs with the same schema as the restricted model inputs.

## Repository modules

1. `data_public/01_city_metadata`: city metadata and contextual indicators.
2. `data_public/02_built_environment_200m`: public 200-m built-environment covariates and functional-context classes.
3. `data_public/03_response_curves`: response-curve exports and nonlinear-model diagnostics.
4. `data_public/04_response_supported_ranges`: response-supported range and threshold sensitivity summaries.
5. `data_public/05_policy_simulation_results`: aggregate policy simulation outputs.

## System requirements

The public figure demo was tested on Windows 10 (22H2) with Python 3.9.7. The
supplied Conda environment uses Python 3.11 and lists its required dependencies.
The full analysis scripts in `code/00_data_preparation` through
`code/07_supplementary` use R; their package requirements are listed in
`environment/R_environment_requirements.md`. The main GAMMs were fitted with `mgcv` 1.9-1,
as reported in the manuscript. No non-standard hardware is required.

## Installation

From the repository root, create and activate the environment:

```bash
conda env create -f environment/environment.yml
conda activate nighttime-pa-lighting-equity-demo
```

Installation typically takes 2-5 minutes on a standard desktop with Conda
already installed.

## Reproducing the figures

Run the release-safe workflow from the repository root `Figure_code/` 
The workflow uses the supplied aggregate tables to reproduce the public
versions of analytic results and Figs. 1-5. Each module writes its PNG and SVG outputs to its own
`Figure_code/` directory. Expected runtime is a few minutes on a standard
desktop.

## Analysis code

The analysis scripts are provided in `code.rar`; `code/00_data_preparation` through `code/07_supplementary`
document the data preparation, descriptive analysis, multilevel models, GAMMs,
response-range extraction, policy simulation and validation workflow. They
require authorized local access to the restricted Strava and SDGSAT-1 inputs
and are therefore not part of the one-command public demo.

The scripts in `Figure_code/` reproduce the public analytic results and figures from the supplied
aggregate and synthetic tables. The descriptive figure modules use city-level
derived summaries; restricted road-level PA records are not redistributed.

## Using other data

Authorized users can prepare a grid-day CSV using the columns and units listed
in `data_demo/restricted_model_input_schema.csv`, then point the corresponding
analysis script to that file. Restricted Strava and SDGSAT-1 data are not
included in this repository.
