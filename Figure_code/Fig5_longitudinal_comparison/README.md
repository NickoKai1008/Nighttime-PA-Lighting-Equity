# Figure 5

Women's nighttime physical activity as daylight fades and lighting changes.

From this directory, run the following commands to write `output/fig5.png` (600 dpi), `output/fig5.svg` and `output/fig5.pdf`.

```sh
python -m pip install numpy==2.3.5 pandas==2.3.3 matplotlib==3.10.6 scipy==1.16.3 Pillow==12.0.0
python plot.py
```

Tested with Python 3.13.9. The script reads the tables and photographs from `data/`. CSV filenames begin with the corresponding panel letter.

| Panel | Source data |
| --- | --- |
| a | Melbourne daily dark minutes, fitted darkness associations for Melbourne and the eight-city pool, the radiance spline and equal-increment contrasts |
| b | Quarterly and annual Hong Kong treatment–control radiance contrasts, paired-window difference-in-differences, and the annotated image pair from On Sau Road Park. The inset labels and dashed lines identify the 2019 pre-renewal and 2023 post-renewal analysis windows. The photographs were acquired in September 2019 and January 2024, respectively (`b_image_pair.csv`). |
| c | Before- and after-renewal streetscape composition across eight Hong Kong sites |
| d | Day and night differences-in-differences for individual projects and equal-project averages |
| e–f | Activity and female-share estimates with 95% confidence intervals and 2,000 matched spatial placebo estimates per contrast |
| g | Sex-specific project ratios, fitted lines and slope estimates for 33 projects across four cities |
| h | Estimated female activity slope by starting radiance, interaction coefficient, and project starting radiances |
