from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CURVE_CSV = DATA_DIR / "policy_curve_profiles.csv"
MAP_CSV = DATA_DIR / "city_action_units.csv"

CITY_ORDER = [
    "Hong Kong",
    "Dubai",
    "Sydney",
    "Auckland",
    "Kuala Lumpur",
    "Doha",
    "Bangkok",
    "Tehran",
    "Dhaka",
]

ACTION_ORDER = [
    "unchanged",
    "dim_noharm",
    "brighten_existing_female_gain",
    "brighten_male_only_to_mixed",
    "brighten_no_activity_to_female_enabled",
]

ACTION_LABELS = {
    "unchanged": "Unchanged / not selected",
    "dim_noharm": "No-harm dimming",
    "brighten_existing_female_gain": "Brightening: existing female PA",
    "brighten_male_only_to_mixed": "Brightening: male-only to mixed",
    "brighten_no_activity_to_female_enabled": "Brightening: no-PA to female-enabled",
}

ACTION_COLORS = {
    "unchanged": "#BEBEBE",
    "dim_noharm": "#8560A8",
    "brighten_existing_female_gain": "#D96B28",
    "brighten_male_only_to_mixed": "#27A7A0",
    "brighten_no_activity_to_female_enabled": "#A9DCD6",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "svg.fonttype": "none",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.dpi": 220,
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
        }
    )


def save_figure(fig: plt.Figure, stem: str, dpi: int = 450) -> None:
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def smooth_series(x: np.ndarray, y: np.ndarray, x_new: np.ndarray, bw: float, weight: np.ndarray) -> np.ndarray:
    weight = np.sqrt(np.maximum(np.asarray(weight, dtype=float), 1.0))
    out = []
    for x0 in x_new:
        w = np.exp(-0.5 * ((x - x0) / bw) ** 2) * weight
        ok = np.isfinite(y) & np.isfinite(w) & (w > 1e-12)
        out.append(np.sum(w[ok] * y[ok]) / np.sum(w[ok]))
    return np.asarray(out)


def x_stretch(x) -> np.ndarray:
    return np.log1p(np.asarray(x, dtype=float) / 0.00055)


def set_curve_ticks(ax: plt.Axes) -> None:
    ticks = np.array([0, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012])
    ax.set_xticks(x_stretch(ticks))
    ax.set_xticklabels([f"{v:0.4f}" if v < 0.002 else f"{v:0.3f}" for v in ticks])


def curve_lines(stats: pd.DataFrame):
    x = stats["x_mid"].to_numpy(dtype=float)
    n = stats["n_grids"].to_numpy(dtype=float)
    x_dense = np.linspace(x.min(), x.max(), 420)
    before_intensity = smooth_series(x, stats["before_intensity"].to_numpy(dtype=float), x_dense, 0.00070, n)
    after_intensity = smooth_series(x, stats["after_intensity"].to_numpy(dtype=float), x_dense, 0.00070, n)
    before_footprint = smooth_series(x, stats["before_footprint"].to_numpy(dtype=float), x_dense, 0.00095, n)
    after_footprint = smooth_series(x, stats["after_footprint"].to_numpy(dtype=float), x_dense, 0.00095, n)
    return x_stretch(x_dense), before_intensity, after_intensity, before_footprint, after_footprint


def plot_policy_curve() -> None:
    stats = pd.read_csv(CURVE_CSV)

    purple = "#B22A8E"
    purple_before = "#8C6A84"
    purple_fill = "#D85AAE"
    teal = "#0E9693"
    teal_before = "#496F73"
    teal_fill = "#4FC4BE"

    x = stats["x_mid"].to_numpy(dtype=float)
    x_plot = x_stretch(x)
    x_dense, before_i, after_i, before_f, after_f = curve_lines(stats)

    fig, ax = plt.subplots(figsize=(7.4, 4.15))
    ax2 = ax.twinx()

    ax.fill_between(x_dense, before_i, after_i, where=after_i >= before_i, color=purple_fill, alpha=0.30, linewidth=0)
    ax2.fill_between(x_dense, before_f, after_f, where=after_f >= before_f, color=teal_fill, alpha=0.24, linewidth=0)

    line1, = ax.plot(x_dense, before_i, color=purple_before, lw=1.8, label="Female PA intensity, before")
    line2, = ax.plot(x_dense, after_i, color=purple, lw=2.55, label="Female PA intensity, after")
    line3, = ax2.plot(x_dense, before_f, color=teal_before, lw=1.75, ls=(0, (4, 3)), label="Female footprint, before")
    line4, = ax2.plot(x_dense, after_f, color=teal, lw=2.35, ls=(0, (4, 3)), label="Female footprint, after")

    selected = stats["bin"].isin([0, 1, 2, 3, 4, 5, 6, 8, 10, 14, 20, 28, 41]).to_numpy()
    for idx in np.where(selected)[0]:
        ax.plot(
            [x_plot[idx], x_plot[idx]],
            [stats["before_intensity"].iloc[idx], stats["after_intensity"].iloc[idx]],
            color=purple,
            lw=1.25,
            alpha=0.78,
            solid_capstyle="round",
            zorder=5,
        )
        ax2.plot(
            [x_plot[idx], x_plot[idx]],
            [stats["before_footprint"].iloc[idx], stats["after_footprint"].iloc[idx]],
            color=teal,
            lw=1.1,
            alpha=0.68,
            solid_capstyle="round",
            zorder=5,
        )

    ax.set_xlim(x_stretch([0])[0], x_stretch([0.0126])[0])
    ax.set_ylim(0.35, 5.05)
    ax2.set_ylim(16, 42)
    set_curve_ticks(ax)
    ax.set_xlabel(r"Baseline ALAN radiance, stretched axis (W m$^{-2}$ sr$^{-1}$ $\mu$m$^{-1}$)")
    ax.set_ylabel("Female PA intensity\n(counts per grid)", color=purple, labelpad=7)
    ax2.set_ylabel("Female-enabled footprint\n(% of grids)", color=teal, labelpad=7)
    ax.tick_params(axis="y", labelcolor=purple)
    ax2.tick_params(axis="y", labelcolor=teal)
    ax.grid(False)
    ax2.grid(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_linewidth(0.8)
    ax.legend(handles=[line1, line2, line3, line4], loc="upper left", bbox_to_anchor=(0.015, 1.02), ncol=2, fontsize=7.1, handlelength=2.5)
    ax.set_title(
        "Reallocation lifts female activity and expands the female-enabled footprint",
        loc="left",
        fontsize=10,
        fontweight="bold",
        pad=15,
    )
    fig.text(
        0.11,
        0.235,
        "Vertical segments and shaded ribbons show the after-before lift on the same curve scale.",
        color="#4B5563",
        fontsize=7.0,
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )
    fig.subplots_adjust(left=0.115, right=0.885, bottom=0.18, top=0.80)
    save_figure(fig, "policy_curve_gain_visible")


def add_local_positions(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["x_km"] = np.nan
    out["y_km"] = np.nan
    for city, idx in out.groupby("display_city").groups.items():
        sub = out.loc[idx]
        lat0 = np.nanmean(sub["lat"])
        lon0 = np.nanmean(sub["lon"])
        out.loc[idx, "x_km"] = (sub["lon"] - lon0) * 111.32 * math.cos(math.radians(lat0))
        out.loc[idx, "y_km"] = (sub["lat"] - lat0) * 110.57
    return out


def unit_polygons(sub: pd.DataFrame, side_km: float = 0.20):
    half = side_km / 2
    return [
        [(x - half, y - half), (x + half, y - half), (x + half, y + half), (x - half, y + half)]
        for x, y in zip(sub["x_km"], sub["y_km"])
    ]


def plot_action_map() -> None:
    data = pd.read_csv(MAP_CSV)
    data["display_city"] = pd.Categorical(data["display_city"], categories=CITY_ORDER, ordered=True)
    data = data.sort_values(["display_city", "lat", "lon"])
    data = add_local_positions(data)

    fig, axes = plt.subplots(3, 3, figsize=(9.0, 8.6))
    axes = axes.flatten()
    for ax, city in zip(axes, CITY_ORDER):
        sub = data[data["display_city"] == city]
        for action in ACTION_ORDER:
            part = sub[sub["action_unified"] == action]
            if len(part) == 0:
                continue
            layer = PolyCollection(
                unit_polygons(part),
                facecolors=ACTION_COLORS[action],
                edgecolors="none",
                linewidths=0,
                alpha=0.96,
                rasterized=True,
                zorder=ACTION_ORDER.index(action),
            )
            ax.add_collection(layer)

        x_min, x_max = sub["x_km"].min(), sub["x_km"].max()
        y_min, y_max = sub["y_km"].min(), sub["y_km"].max()
        span = max(x_max - x_min, y_max - y_min, 1.0)
        pad = max(span * 0.035, 0.3)
        x_mid = (x_min + x_max) / 2
        y_mid = (y_min + y_max) / 2
        ax.set_xlim(x_mid - span / 2 - pad, x_mid + span / 2 + pad)
        ax.set_ylim(y_mid - span / 2 - pad, y_mid + span / 2 + pad)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(city, fontsize=8, pad=3)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    handles = [Patch(facecolor=ACTION_COLORS[k], edgecolor="none", label=ACTION_LABELS[k]) for k in ACTION_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle("Policy action map across selected cities", fontsize=12, fontweight="bold", y=0.995)
    fig.subplots_adjust(left=0.025, right=0.985, top=0.95, bottom=0.12, wspace=0.04, hspace=0.13)
    save_figure(fig, "selected_city_action_change_map", dpi=360)


def main() -> None:
    setup_style()
    plot_policy_curve()
    plot_action_map()
    print("Saved figures to", FIG_DIR)


if __name__ == "__main__":
    main()
