from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CURVE_CSV = DATA_DIR / "lorenz_curve_points.csv"
GINI_CSV = DATA_DIR / "night_gini_by_city.csv"

GII_ORDER = ["Low", "Medium", "High"]
REGION_ORDER = [
    "East Asia",
    "Southeast Asia",
    "South & Central Asia",
    "Pacific",
    "Middle East",
]

GII_BASE = {
    "Low": "#FF8C2A",
    "Medium": "#27B7AE",
    "High": "#8D6AD8",
}
REGION_COLORS = {
    "East Asia": "#3C5488",
    "Southeast Asia": "#00A087",
    "South & Central Asia": "#E64B35",
    "Pacific": "#4DBBD5",
    "Middle East": "#F39B7F",
}


def setup_style():
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "svg.fonttype": "none",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def lighten(hex_color, amount=0.55):
    hex_color = hex_color.lstrip("#")
    rgb = np.array([int(hex_color[i:i + 2], 16) for i in (0, 2, 4)]) / 255
    out = rgb + (1 - rgb) * amount
    return tuple(out)


def group_gradient_color(group, gini, group_min, group_max):
    base = np.array(mpl.colors.to_rgb(GII_BASE[group]))
    if group_max <= group_min:
        weight = 0.55
    else:
        weight = (gini - group_min) / (group_max - group_min)
    pale = np.array(lighten(GII_BASE[group], 0.62))
    return tuple(pale * (1 - weight) + base * weight)


def plot_lorenz(curves, outcome, stem, ylabel):
    subset = curves[curves["outcome"].eq(outcome)].copy()
    city_meta = (
        subset[["city", "city_key", "gii_group", "gini"]]
        .drop_duplicates()
        .sort_values("gini")
    )

    ranges = city_meta.groupby("gii_group")["gini"].agg(["min", "max"]).to_dict("index")

    fig = plt.figure(figsize=(18.0, 8.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.42, 0.88], wspace=0.22)
    ax_main = fig.add_subplot(gs[0, 0])
    ax_legend = fig.add_subplot(gs[0, 1])
    ax_zoom = fig.add_subplot(gs[0, 2])
    ax_legend.axis("off")

    ax_main.plot([0, 1], [0, 1], ls="--", color="#777777", lw=1.2, label="Perfect equality", zorder=1)
    ax_zoom.plot([0.5, 1], [0, 0.5], ls="--", color="#777777", lw=1.2, zorder=1)

    for _, meta in city_meta.iterrows():
        city_curve = subset[subset["city_key"].eq(meta["city_key"])]
        group = meta["gii_group"]
        color = group_gradient_color(group, meta["gini"], ranges[group]["min"], ranges[group]["max"])
        label = f"{meta['city']} ({meta['gini']:.3f}) [{group[0]}]"
        ax_main.plot(
            city_curve["population_fraction"],
            city_curve["activity_fraction"],
            color=color,
            lw=1.35,
            alpha=0.88,
            label=label,
            zorder=2,
        )
        ax_zoom.plot(
            city_curve["population_fraction"],
            city_curve["activity_fraction"],
            color=color,
            lw=1.5,
            alpha=0.88,
            zorder=2,
        )

    ax_main.set_xlim(0, 1)
    ax_main.set_ylim(0, 1)
    ax_main.set_aspect("equal", adjustable="box")
    ax_main.set_xlabel("Cumulative fraction of spatial units", fontsize=10.5, fontweight="bold")
    ax_main.set_ylabel(ylabel, fontsize=10.5, fontweight="bold")
    ax_main.set_title(ylabel.replace("Cumulative fraction of ", ""), fontsize=12, fontweight="bold", pad=10)
    ax_main.tick_params(labelsize=8.5, length=3)

    ax_zoom.set_xlim(0.5, 1.0)
    ax_zoom.set_ylim(0.0, 0.5)
    ax_zoom.set_aspect("equal", adjustable="box")
    ax_zoom.set_xlabel("Cumulative fraction, zoomed", fontsize=10, fontweight="bold")
    ax_zoom.set_ylabel("Cumulative fraction, zoomed", fontsize=10, fontweight="bold")
    ax_zoom.set_title("Right-bottom corner zoom", fontsize=11, fontweight="bold", pad=10)
    ax_zoom.tick_params(labelsize=8.5, length=3)

    counts = city_meta["gii_group"].value_counts().reindex(GII_ORDER).fillna(0).astype(int)
    note = "GII classification:\n" + "\n".join([f"{g}: {counts[g]} cities" for g in GII_ORDER])
    ax_main.text(
        0.03, 0.35, note, fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#333333", linewidth=0.8),
    )

    handles, labels = ax_main.get_legend_handles_labels()
    ax_legend.legend(
        handles, labels, loc="center", fontsize=5.7,
        title="City (Gini)", title_fontsize=8, frameon=True,
        edgecolor="#333333", facecolor="white", fancybox=False,
    )

    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.12, top=0.90)
    for ext in ["png", "svg"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def box_kwargs():
    return dict(
        patch_artist=True,
        showfliers=False,
        showmeans=True,
        meanprops=dict(marker="D", mfc="white", mec="#333333", ms=4.3, zorder=4),
        medianprops=dict(color="#111111", lw=2.0, solid_capstyle="round"),
        whiskerprops=dict(color="#555555", lw=1.0),
        capprops=dict(color="#555555", lw=1.0),
    )


def plot_boxplot(gini):
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(7.5, 4.0), sharey=True,
        gridspec_kw={"width_ratios": [3, 5], "wspace": 0.05},
    )
    y = gini["night_female_gini"]
    pad = (y.max() - y.min()) * 0.12
    y_lo = max(0, y.min() - pad)
    y_hi = min(1.02, y.max() + pad)

    gii_data = [gini.loc[gini["gii_group"].eq(g), "night_female_gini"].values for g in GII_ORDER]
    bp_a = ax_a.boxplot(gii_data, positions=range(len(GII_ORDER)), widths=0.58, **box_kwargs())
    for box, group in zip(bp_a["boxes"], GII_ORDER):
        box.set_facecolor(lighten(GII_BASE[group], 0.62))
        box.set_edgecolor(GII_BASE[group])
        box.set_linewidth(1.6)

    ax_a.set_xticks(range(len(GII_ORDER)))
    ax_a.set_xticklabels([f"{g}\n(n = {len(gii_data[i])})" for i, g in enumerate(GII_ORDER)])
    ax_a.set_xlabel("UNDP Gender Inequality Index", fontsize=9, labelpad=8)
    ax_a.set_ylabel("Gini coefficient", fontsize=9, labelpad=8)
    ax_a.set_ylim(y_lo, y_hi)
    ax_a.text(-0.78, y_hi + (y_hi - y_lo) * 0.04, "a", fontsize=12, fontweight="bold", va="bottom", clip_on=False)

    region_data = [gini.loc[gini["region"].eq(r), "night_female_gini"].values for r in REGION_ORDER]
    bp_b = ax_b.boxplot(region_data, positions=range(len(REGION_ORDER)), widths=0.58, **box_kwargs())
    for box, region in zip(bp_b["boxes"], REGION_ORDER):
        box.set_facecolor(lighten(REGION_COLORS[region], 0.62))
        box.set_edgecolor(REGION_COLORS[region])
        box.set_linewidth(1.6)

    short = ["East\nAsia", "Southeast\nAsia", "South &\nCentral Asia", "Pacific", "Middle\nEast"]
    ax_b.set_xticks(range(len(REGION_ORDER)))
    ax_b.set_xticklabels([f"{s}\n(n = {len(region_data[i])})" for i, s in enumerate(short)])
    ax_b.set_xlabel("Geographic region", fontsize=9, labelpad=8)
    ax_b.text(-0.78, y_hi + (y_hi - y_lo) * 0.04, "b", fontsize=12, fontweight="bold", va="bottom", clip_on=False)

    for ax in (ax_a, ax_b):
        ax.tick_params(axis="x", length=0, pad=4)
        ax.tick_params(axis="y", labelsize=8)
        ax.yaxis.grid(True, ls="-", lw=0.35, alpha=0.25, color="#AAAAAA")
        ax.set_axisbelow(True)

    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.24, top=0.90, wspace=0.05)
    for ext in ["png", "svg"]:
        fig.savefig(FIG_DIR / f"night_female_gini_boxplot.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    setup_style()
    curves = pd.read_csv(CURVE_CSV)
    gini = pd.read_csv(GINI_CSV)
    plot_lorenz(
        curves,
        outcome="night_female",
        stem="lorenz_night_female",
        ylabel="Cumulative fraction of nighttime female activity",
    )
    plot_lorenz(
        curves,
        outcome="night_male",
        stem="lorenz_night_male",
        ylabel="Cumulative fraction of nighttime male activity",
    )
    plot_boxplot(gini)
    manifest = "\n".join(str(p.relative_to(ROOT)).replace("\\", "/") for p in sorted(FIG_DIR.glob("*")))
    (ROOT / "logs" / "figure_manifest.txt").write_text(manifest, encoding="utf-8")


if __name__ == "__main__":
    main()
