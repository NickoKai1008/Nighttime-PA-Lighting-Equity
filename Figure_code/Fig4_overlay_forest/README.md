# Right-label overlay forest figures

This folder contains a clean, reproducible package for the final female-male overlay forest figures with right-aligned coefficient and confidence-interval labels.

## Files

- `scripts/overlay_forest.py`: standalone plotting script.
- `data/lmm_plot_data.csv`: figure-ready LMM coefficient table for female and male nighttime PA.
- `data/gamm_plot_data.csv`: figure-ready adjusted GAMM parametric coefficient table for female and male nighttime PA.
- `raw_model_results/`: source model-result CSVs used to build the figure-ready tables.
- `figures/`: regenerated PNG and editable SVG outputs.

## Run

```powershell
.\Fig4_overlay_forest\scripts\overlay_forest.py
```

The script reads only local files from `data/` and writes to `figures/`.
