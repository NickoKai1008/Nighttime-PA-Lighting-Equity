from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"

PURPLE = "#8F3FA3"
PURPLE_FILL = "#DCC7E1"
TEAL = "#159CA2"
TEAL_LIGHT = "#DDF2F0"
GREEN = "#35A96B"
BLUE = "#3D82C4"
ORANGE = "#F0782A"
GREY = "#DADADA"
BLACK = "#222222"

LANDUSE_ORDER = ["recreational", "residential", "commercial"]
LANDUSE_LABEL = {
    "recreational": "Recreational",
    "residential": "Residential",
    "commercial": "Commercial",
}
LANDUSE_COLOR = {
    "recreational": GREEN,
    "residential": BLUE,
    "commercial": ORANGE,
}

def set_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "svg.fonttype": "none",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.1,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "legend.frameon": False,
            "figure.dpi": 160,
        }
    )


def save_figure(fig, stem):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.png", dpi=420, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def robust_ylim(y, lo, hi, pad=0.08):
    vals = np.concatenate([np.asarray(y), np.asarray(lo), np.asarray(hi)])
    vals = vals[np.isfinite(vals)]
    y0 = np.nanpercentile(vals, 1)
    y1 = np.nanpercentile(vals, 99)
    span = max(y1 - y0, 1e-6)
    return y0 - pad * span, y1 + pad * span


def robust_xlim(*arrays, pad=0.03):
    vals = np.concatenate([np.asarray(a, dtype=float).ravel() for a in arrays if len(a)])
    vals = vals[np.isfinite(vals)]
    x0 = float(np.nanmin(vals))
    x1 = float(np.nanmax(vals))
    span = max(x1 - x0, 1e-9)
    return x0 - pad * span, x1 + pad * span


def draw_scaled_hist(ax, hist, hist_base, hist_height, color=GREY):
    if hist.empty:
        return
    h = hist.sort_values("bin_mid")
    max_density = float(h["density"].max())
    if max_density <= 0:
        return
    ax.bar(
        h["bin_left"].to_numpy(float),
        h["density"].to_numpy(float) / max_density * hist_height,
        width=(h["bin_right"] - h["bin_left"]).to_numpy(float),
        bottom=hist_base,
        align="edge",
        color=color,
        edgecolor="#A7A7A7",
        linewidth=0.25,
        zorder=0,
    )


def draw_stacked_hist(ax, hist, hist_base, hist_height):
    h0 = hist[hist["landuse"].notna()].copy()
    if h0.empty:
        return
    totals = h0.groupby("bin_left", observed=True)["density"].sum()
    max_total = float(totals.max()) if len(totals) else 0.0
    if max_total <= 0:
        return
    bottoms = {}
    for lu in LANDUSE_ORDER:
        hh = h0[h0["landuse"].eq(lu)].sort_values("bin_left")
        for _, r in hh.iterrows():
            x0 = float(r["bin_left"])
            width = float(r["bin_right"] - r["bin_left"])
            height = float(r["density"]) / max_total * hist_height
            bottom = bottoms.get(x0, hist_base)
            ax.bar(
                x0,
                height,
                width=width,
                align="edge",
                bottom=bottom,
                color=LANDUSE_COLOR[lu],
                edgecolor="none",
                alpha=0.25,
                zorder=0,
            )
            bottoms[x0] = bottom + height


def add_thresholds(ax, low, peak, color=PURPLE, shade=True, lw=1.35):
    if shade and np.isfinite(low) and np.isfinite(peak):
        ax.axvspan(min(low, peak), max(low, peak), color=TEAL_LIGHT, alpha=0.65, lw=0, zorder=-1)
    if np.isfinite(low):
        ax.axvline(low, color=TEAL, lw=1.05, ls=(0, (2, 2)), zorder=3)
    if np.isfinite(peak):
        ax.axvline(peak, color=color, lw=lw, zorder=3)


def panel_a(ax):
    response = pd.read_csv(DATA / "panel_a_pooled_response.csv").sort_values("ntl_raw")
    thresholds = pd.read_csv(DATA / "panel_a_pooled_thresholds.csv").iloc[0]
    hist = pd.read_csv(DATA / "observed_ntl_histograms_female_dropzero_log_raw.csv")
    hist = hist[hist["landuse"].isna()]

    x = response["ntl_raw"].to_numpy(float)
    y = response["fit"].to_numpy(float)
    lo = response["lo"].to_numpy(float)
    hi = response["hi"].to_numpy(float)
    y0, y1 = robust_ylim(y, lo, hi, pad=0.07)
    span = y1 - y0
    hist_base = y0 - 0.26 * span
    hist_height = 0.17 * span

    low = float(thresholds["operating_low_raw"])
    peak = float(thresholds["derivative_peak_raw"])

    hist_edges = pd.concat([hist["bin_left"], hist["bin_right"]]).to_numpy(float)
    xmin, xmax = robust_xlim(x, hist_edges, [low, peak])

    ax.axvspan(peak, xmax, color="#E8E8E8", alpha=0.88, zorder=-4)
    ax.axhline(y0, color="#9C9C9C", lw=0.65, zorder=1)
    ax.fill_between(x, lo, hi, color=PURPLE_FILL, alpha=0.75, lw=0)
    ax.plot(x, y, color=PURPLE, lw=2.0)
    add_thresholds(ax, low, peak, color=PURPLE, shade=False, lw=1.65)
    draw_scaled_hist(ax, hist, hist_base, hist_height, color="#D9D9D9")

    ax.scatter([peak], [float(np.interp(peak, x, y))], s=22, color=PURPLE, zorder=4)
    ax.text(
        0.58,
        0.50,
        "Energy saving\n&\nlight pollution\nmitigation\npotential",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.2,
        color=BLACK,
    )
    ax.text(0.96, 0.035, "NTL distribution", transform=ax.transAxes, ha="right", va="center", fontsize=7, color="#8A8A8A")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(hist_base, y1)
    ax.set_title("Female nighttime PA ~ NTL (all cities)", fontsize=10.5, y=1.14, pad=0)
    ax.text(0.5, 1.035, "Female nighttime PA", transform=ax.transAxes, ha="center", va="bottom", fontsize=11.5, weight="bold")
    ax.set_ylabel("Female Nighttime PA (log)", fontsize=10.5)
    ax.set_xlabel("NTL radiance (Wm$^{-2}$sr$^{-1}$um$^{-1}$)", fontsize=9.8)
    ax.tick_params(labelsize=8.5)


def panel_c(ax):
    response = pd.read_csv(DATA / "panel_c_landuse_response.csv")
    thresholds = pd.read_csv(DATA / "panel_c_landuse_thresholds.csv").set_index("landuse")
    hist = pd.read_csv(DATA / "observed_ntl_histograms_female_dropzero_log_raw.csv")

    vals = response["fit"].to_numpy(float)
    lows = response["lo"].to_numpy(float)
    highs = response["hi"].to_numpy(float)
    y0, y1 = robust_ylim(vals, lows, highs, pad=0.07)
    span = y1 - y0
    hist_base = y0 - 0.25 * span
    hist_height = 0.18 * span
    ax.axhline(y0, color="#9C9C9C", lw=0.65, zorder=1)

    handles = []
    for lu in LANDUSE_ORDER:
        d = response[response["landuse"].eq(lu)].sort_values("ntl_raw")
        th = thresholds.loc[lu]
        color = LANDUSE_COLOR[lu]
        x = d["ntl_raw"].to_numpy(float)
        y = d["fit"].to_numpy(float)
        lo = d["lo"].to_numpy(float)
        hi = d["hi"].to_numpy(float)
        ax.fill_between(x, lo, hi, color=color, alpha=0.12, lw=0)
        ax.plot(x, y, color=color, lw=1.7)
        add_thresholds(ax, float(th["operating_low_raw"]), float(th["derivative_peak_raw"]), color=color, shade=False, lw=1.2)
        handles.append(Line2D([0], [0], color=color, lw=1.8, label=LANDUSE_LABEL[lu]))

    hist_edges = pd.concat([hist["bin_left"], hist["bin_right"]]).to_numpy(float)
    thresh_vals = thresholds[["operating_low_raw", "derivative_peak_raw"]].to_numpy(float).ravel()
    xmin, xmax = robust_xlim(response["ntl_raw"].to_numpy(float), hist_edges, thresh_vals)

    draw_stacked_hist(ax, hist, hist_base, hist_height)
    ax.text(0.98, 0.035, "stacked NTL distribution", transform=ax.transAxes, ha="right", va="center", fontsize=7, color="#8A8A8A")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(hist_base, y1)
    ax.set_title("Female nighttime PA ~ NTL (by land-use)", fontsize=10.5, pad=10)
    ax.set_ylabel("Adjusted partial effect", fontsize=9.4)
    ax.set_xlabel("NTL radiance (Wm$^{-2}$sr$^{-1}$um$^{-1}$)", fontsize=9.3)
    ax.legend(handles=handles, ncol=3, loc="upper right", fontsize=7.1, handlelength=1.8, columnspacing=1.3)
    ax.tick_params(labelsize=8.2)


def panel_b_grid(fig, spec):
    response = pd.read_csv(DATA / "panel_b_city_conditioned_response_9cities.csv")
    thresholds = pd.read_csv(DATA / "panel_b_city_thresholds_9cities.csv").set_index("city")
    city_order = (
        response[["city", "city_order"]]
        .drop_duplicates()
        .sort_values("city_order")["city"]
        .tolist()
    )
    gs = spec.subgridspec(3, 3, wspace=0.25, hspace=0.45)
    axes = []
    for idx, city in enumerate(city_order):
        ax = fig.add_subplot(gs[idx // 3, idx % 3])
        axes.append(ax)
        d = response[response["city"].eq(city)].sort_values("ntl_raw")
        th = thresholds.loc[city]
        x = d["ntl_raw"].to_numpy(float)
        y = d["fit"].to_numpy(float)
        lo = d["lo"].to_numpy(float)
        hi = d["hi"].to_numpy(float)
        low = float(th["operating_low_raw"])
        peak = float(th["derivative_peak_raw"])

        ax.fill_between(x, lo, hi, color=PURPLE_FILL, alpha=0.58, lw=0)
        ax.plot(x, y, color=PURPLE, lw=1.45)
        add_thresholds(ax, low, peak, color=PURPLE, shade=True, lw=1.15)
        ax.set_title(str(th["city_label"]), fontsize=8, weight="bold", pad=3)
        ax.set_ylim(*robust_ylim(y, lo, hi, pad=0.09))
        xmin, xmax = float(np.nanmin(x)), float(np.nanmax(x))
        ax.set_xlim(xmin, xmax)
        ax.tick_params(labelsize=6.7, length=2.2)
        if idx < 6:
            ax.set_xticklabels([])
    return axes


def make_composite():
    set_style()
    fig = plt.figure(figsize=(13.2, 14.6))
    outer = GridSpec(2, 2, figure=fig, height_ratios=[0.95, 1.86], width_ratios=[1.0, 1.35], hspace=0.30, wspace=0.18)
    ax_a = fig.add_subplot(outer[0, 0])
    ax_c = fig.add_subplot(outer[0, 1])
    panel_a(ax_a)
    panel_c(ax_c)
    axes_b = panel_b_grid(fig, outer[1, :])

    fig.text(0.020, 0.965, "a", fontsize=25, weight="bold")
    fig.text(0.485, 0.965, "c", fontsize=25, weight="bold")
    fig.text(0.020, 0.545, "b", fontsize=25, weight="bold")
    fig.text(0.018, 0.305, "Female Nighttime PA (log)", rotation=90, va="center", ha="center", fontsize=11)
    fig.text(0.505, 0.030, "Nighttime light radiance", va="center", ha="center", fontsize=11)
    save_figure(fig, "fig3_nonlinear_response_composite")


def make_standalone_panels():
    set_style()
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    panel_a(ax)
    save_figure(fig, "fig3_panel_a_pooled_response")

    fig = plt.figure(figsize=(10.2, 8.6))
    gs = GridSpec(1, 1, figure=fig)
    panel_b_grid(fig, gs[0])
    fig.text(0.02, 0.50, "Female Nighttime PA (log)", rotation=90, va="center", ha="center", fontsize=11)
    fig.text(0.50, 0.02, "Nighttime light radiance", va="center", ha="center", fontsize=11)
    save_figure(fig, "fig3_panel_b_city_conditioned_9cities")

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    panel_c(ax)
    save_figure(fig, "fig3_panel_c_landuse_overlay")


def main():
    make_composite()
    make_standalone_panels()


if __name__ == "__main__":
    main()
