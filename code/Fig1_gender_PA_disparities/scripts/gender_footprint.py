from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.stats import t as student_t
except Exception:  # pragma: no cover
    gaussian_filter1d = None
    student_t = None


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
MAP_DIR = ROOT / "maps"
LOG_DIR = ROOT / "logs"
for folder in [FIG_DIR, MAP_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

FEMALE = "#7B3F98"
MALE = "#E6863B"
MIXED = "#1B9E8A"
EMPTY = "#D8D8D8"
TEXT = "#222222"
CI_GREY = "#BDBDBD"

CITY_LABEL = {
    "HongKong": "Hong Kong",
    "Dhaka": "Dhaka",
    "BandarSeriBegawan": "Bandar Seri Begawan",
    "HoChiMinh": "Ho Chi Minh",
    "KualaLumpur": "Kuala Lumpur",
    "Toshkent": "Tashkent",
}

GII_COLORS = {
    "Low GII": "#E6863B",
    "Medium GII": "#55BFC0",
    "High GII": "#8E63B0",
}

REGION_MARKERS = {
    "East Asia": "o",
    "Southeast Asia": "^",
    "South/Central Asia": "s",
    "Pacific": "X",
    "Middle East": "D",
}

CATEGORY_STYLE = {
    "both_zero": ("No female or male PA", EMPTY, 0.48, 1),
    "both_active": ("Gender-mixed PA", MIXED, 0.90, 2),
    "female_only": ("Female-only PA", FEMALE, 0.96, 3),
    "male_only": ("Male-only PA", MALE, 0.92, 4),
}


def apply_style(font_size: float = 8) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "svg.fonttype": "none",
            "font.size": font_size,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
        }
    )


def city_display(city: str) -> str:
    return CITY_LABEL.get(city, city)


def add_scale_bar(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray, length_km: float = 5.0) -> None:
    if len(lon) == 0 or len(lat) == 0:
        return
    lat0 = float(np.nanmedian(lat))
    km_per_degree_lon = 111.32 * max(math.cos(math.radians(lat0)), 0.15)
    bar_deg = length_km / km_per_degree_lon
    x0 = float(np.nanmin(lon) + 0.07 * (np.nanmax(lon) - np.nanmin(lon)))
    y0 = float(np.nanmin(lat) + 0.07 * (np.nanmax(lat) - np.nanmin(lat)))
    ax.plot([x0, x0 + bar_deg], [y0, y0], color="#333333", lw=1.0, solid_capstyle="butt")
    ax.text(x0, y0 + 0.012 * (np.nanmax(lat) - np.nanmin(lat)), f"{int(length_km)} km", ha="left", va="bottom", fontsize=6.5, color=TEXT)


def plot_maps(grid: pd.DataFrame) -> None:
    apply_style(font_size=8)
    cities = ["HongKong", "Dhaka"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), constrained_layout=True)
    for ax, city in zip(axes, cities):
        sub = grid[grid["city"].eq(city)].copy()
        ax.set_title(city_display(city), loc="left", fontsize=10, fontweight="bold", pad=2)
        for category, (label, color, alpha, zorder) in CATEGORY_STYLE.items():
            dd = sub[sub["zero_category"].eq(category)]
            if dd.empty:
                continue
            size = 3.8 if city == "HongKong" else 5.0
            ax.scatter(dd["lon"], dd["lat"], s=size, marker="s", c=color, alpha=alpha, edgecolors="none", rasterized=False, zorder=zorder)
        lon = sub["lon"].to_numpy(dtype=float)
        lat = sub["lat"].to_numpy(dtype=float)
        pad_x = max((np.nanmax(lon) - np.nanmin(lon)) * 0.035, 0.005)
        pad_y = max((np.nanmax(lat) - np.nanmin(lat)) * 0.035, 0.005)
        ax.set_xlim(np.nanmin(lon) - pad_x, np.nanmax(lon) + pad_x)
        ax.set_ylim(np.nanmin(lat) - pad_y, np.nanmax(lat) + pad_y)
        ax.set_aspect(1 / max(math.cos(math.radians(float(np.nanmedian(lat)))), 0.15))
        add_scale_bar(ax, lon, lat, length_km=5.0)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("#D0D0D0")
            spine.set_linewidth(0.55)
    handles = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=color, markeredgecolor="none", markersize=7, label=label, alpha=alpha)
        for _, (label, color, alpha, _) in CATEGORY_STYLE.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, bbox_to_anchor=(0.5, -0.035), handletextpad=0.35, columnspacing=0.95)
    for ext in ["png", "svg"]:
        fig.savefig(MAP_DIR / f"map_hongkong_dhaka_nighttime_gender_footprint_200m.{ext}", dpi=520 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def smooth_density(values: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals = values[np.isfinite(values)]
    vals = vals[vals >= 0]
    if len(vals) < 3 or np.nanmax(vals) <= np.nanmin(vals):
        return np.array([]), np.array([])
    counts, edges = np.histogram(vals, bins=bins, density=True)
    if gaussian_filter1d is not None:
        counts = gaussian_filter1d(counts.astype(float), sigma=1.15)
    centers = (edges[:-1] + edges[1:]) / 2.0
    return centers, counts


def plot_density(daily: pd.DataFrame) -> None:
    apply_style(font_size=8)
    cities = ["HongKong", "Dhaka"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), constrained_layout=True)
    for ax, city in zip(axes, cities):
        sub = daily[daily["city"].eq(city)]
        female = sub[sub["sex"].eq("female")]["footprint_area_km2"].to_numpy(dtype=float)
        male = sub[sub["sex"].eq("male")]["footprint_area_km2"].to_numpy(dtype=float)
        vals = np.concatenate([female, male])
        upper = max(float(np.nanpercentile(vals, 99.0)), float(np.nanmax(vals)), 1.0)
        bins = np.linspace(0, upper, 24)
        ax.hist(male, bins=bins, density=True, color=MALE, alpha=0.12, edgecolor="none")
        ax.hist(female, bins=bins, density=True, color=FEMALE, alpha=0.18, edgecolor="none")
        for arr, color, lw, alpha_fill in [(male, MALE, 1.05, 0.14), (female, FEMALE, 1.75, 0.25)]:
            x, y = smooth_density(arr, bins)
            if len(x):
                ax.fill_between(x, y, color=color, alpha=alpha_fill, linewidth=0)
                ax.plot(x, y, color=color, lw=lw, alpha=0.9)
        ax.set_xlim(0, upper * 1.02)
        ax.set_title(city_display(city), loc="left", fontsize=10, fontweight="bold", pad=2)
        ax.set_xlabel("Daily nighttime footprint area (km2)", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.grid(axis="y", color="#EAEAEA", linewidth=0.5)
        ax.tick_params(labelsize=7, length=2)
    handles = [
        Line2D([0], [0], color=FEMALE, lw=1.8, label="Female PA footprint"),
        Line2D([0], [0], color=MALE, lw=1.2, label="Male PA footprint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, -0.055), handlelength=2.2)
    for ext in ["png", "svg"]:
        fig.savefig(FIG_DIR / f"density_hongkong_dhaka_daily_nighttime_footprint_area_by_sex.{ext}", dpi=450 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fit_line_with_ci(x: np.ndarray, y: np.ndarray, xx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    slope, intercept = np.polyfit(x, y, 1)
    yy = slope * xx + intercept
    yhat = slope * x + intercept
    residual = y - yhat
    n = len(x)
    s_err = math.sqrt(float(np.sum(residual**2)) / max(n - 2, 1))
    xbar = float(np.mean(x))
    sxx = float(np.sum((x - xbar) ** 2))
    if sxx <= 0:
        se = np.zeros_like(xx)
    else:
        se = s_err * np.sqrt((1.0 / n) + ((xx - xbar) ** 2 / sxx))
    crit = float(student_t.ppf(0.975, n - 2)) if student_t is not None and n > 2 else 2.064
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return yy, yy - crit * se, yy + crit * se, r2


def plot_footprint_share_scatter(scatter: pd.DataFrame) -> None:
    apply_style(font_size=8)
    xcol = "female_footprint_share_night_pct"
    ycol = "female_footprint_share_day_pct"
    d = scatter.dropna(subset=[xcol, ycol]).copy()
    fig, ax = plt.subplots(figsize=(6.25, 5.65))
    lim_max = max(45.0, math.ceil(max(float(d[xcol].max()), float(d[ycol].max())) / 5.0) * 5.0)
    ax.fill_between([0, lim_max], [0, lim_max], [lim_max, lim_max], color="#F0F0F0", alpha=0.72, zorder=0)
    ax.plot([0, lim_max], [0, lim_max], color="#8A8A8A", lw=1.1, ls=(0, (4, 3)), zorder=1)

    x = d[xcol].to_numpy(dtype=float)
    y = d[ycol].to_numpy(dtype=float)
    xx = np.linspace(0, lim_max, 200)
    yy, lo, hi, r2 = fit_line_with_ci(x, y, xx)
    ax.fill_between(xx, lo, hi, color=CI_GREY, alpha=0.28, linewidth=0, zorder=1)
    ax.plot(xx, yy, color=FEMALE, lw=1.55, zorder=2, label=f"Linear fit (R2 = {r2:.2f})")

    for region, marker in REGION_MARKERS.items():
        region_sub = d[d["region"].eq(region)]
        if region_sub.empty:
            continue
        for group, color in GII_COLORS.items():
            sub = region_sub[region_sub["gii_group"].eq(group)]
            if sub.empty:
                continue
            ax.scatter(sub[xcol], sub[ycol], s=50, marker=marker, color=color, alpha=0.94, edgecolor="white", linewidth=0.65, zorder=4)

    for _, row in d.iterrows():
        ax.annotate(city_display(row["city"]), xy=(row[xcol], row[ycol]), xytext=(5, 4), textcoords="offset points", fontsize=5.8, color=TEXT, zorder=5)

    ax.set_xlim(0, lim_max)
    ax.set_ylim(0, lim_max)
    ticks = np.arange(0, lim_max + 0.1, 10.0)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel("Nighttime female footprint share (%)", fontsize=8.7)
    ax.set_ylabel("Daytime female footprint share (%)", fontsize=8.7)
    ax.set_title("Female share of PA footprint area", loc="left", fontsize=10.3, fontweight="bold")
    ax.text(lim_max - 0.45, lim_max - 0.65, "Higher female share during daytime", ha="right", va="top", fontsize=6.6, fontstyle="italic", color="#333333")
    ax.tick_params(axis="both", direction="out", length=4.0, width=0.9, color=TEXT)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(1.0)
        ax.spines[spine].set_color(TEXT)

    group_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="white", markeredgewidth=0.65, markersize=5.8, label=group.replace(" GII", ""))
        for group, color in GII_COLORS.items()
    ]
    region_handles = [
        Line2D([0], [0], marker=marker, color="#666666", linestyle="none", markerfacecolor="#BFBFBF", markeredgecolor="white", markeredgewidth=0.65, markersize=5.8, label=region)
        for region, marker in REGION_MARKERS.items()
    ]
    fit_handle = Line2D([0], [0], color=FEMALE, lw=1.55, label=f"Linear fit (R2 = {r2:.2f})")
    band_handle = Line2D([0], [0], color=CI_GREY, lw=6.0, alpha=0.35, label="95% confidence band")
    diag_handle = Line2D([0], [0], color="#8A8A8A", lw=1.1, ls=(0, (4, 3)), label="Daytime = nighttime")
    leg1 = ax.legend(handles=group_handles, title="UNDP GII", loc="lower right", bbox_to_anchor=(0.995, 0.18), fontsize=6.5, title_fontsize=6.8, ncol=3, handletextpad=0.3, columnspacing=0.55)
    ax.add_artist(leg1)
    ax.legend(handles=region_handles + [fit_handle, band_handle, diag_handle], title="Geographic region", loc="upper left", bbox_to_anchor=(0.02, 0.985), fontsize=6.3, title_fontsize=6.7, handletextpad=0.45, labelspacing=0.35)

    fig.tight_layout()
    stem = "scatter_daytime_vs_nighttime_female_footprint_share_25cities_fig1b_style"
    for ext in ["png", "svg"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=520 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    summary = pd.DataFrame(
        [
            {
                "n_cities": len(d),
                "night_lower_than_day_n": int((d[xcol] < d[ycol]).sum()),
                "median_day": float(np.median(y)),
                "median_night": float(np.median(x)),
                "mean_day": float(np.mean(y)),
                "mean_night": float(np.mean(x)),
                "linear_r2": float(r2),
            }
        ]
    )
    summary.to_csv(DATA_DIR / "female_footprint_share_scatter_summary.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    grid = pd.read_csv(DATA_DIR / "gender_grid_month_categories_200m.csv")
    daily = pd.read_csv(DATA_DIR / "daily_nighttime_footprint_area_by_city_sex.csv")
    scatter = pd.read_csv(DATA_DIR / "city_day_night_female_footprint_share_25cities.csv")
    plot_maps(grid)
    plot_density(daily)
    plot_footprint_share_scatter(scatter)
    files = sorted(list(FIG_DIR.glob("*")) + list(MAP_DIR.glob("*")))
    manifest = "\n".join(path.relative_to(ROOT).as_posix() for path in files)
    (LOG_DIR / "figure_manifest.txt").write_text(manifest, encoding="utf-8")
    print(manifest)


if __name__ == "__main__":
    main()
