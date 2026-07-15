from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "gender_inequality_quadrant_data.csv"
FIG_DIR = ROOT / "figures"
LOG_DIR = ROOT / "logs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "Low GII": "#E6863B",
    "Medium GII": "#55BFC0",
    "High GII": "#8E63B0",
}
MARKERS = {
    "East Asia": "o",
    "Southeast Asia": "^",
    "South/Central Asia": "s",
    "Pacific": "X",
    "Middle East": "D",
}

def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "svg.fonttype": "none",
            "font.size": 11,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 1.0,
            "legend.frameon": False,
        }
    )


def _overlap_area(a, b) -> float:
    width = max(0.0, min(a.x1, b.x1) - max(a.x0, b.x0))
    height = max(0.0, min(a.y1, b.y1) - max(a.y0, b.y0))
    return width * height


def annotate_cities(ax, data: pd.DataFrame, x_col: str, y_col: str) -> None:
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)
    radii = [7, 11, 16, 22, 29, 37, 46]
    angles = np.deg2rad([30, -30, 150, -150, 0, 90, -90, 180])
    candidates = [(r * np.cos(a), r * np.sin(a)) for r in radii for a in angles]

    points = ax.transData.transform(data[[x_col, y_col]].to_numpy(dtype=float))
    scale = np.maximum(np.ptp(points, axis=0), 1.0)
    normalized = points / scale
    distance = np.sqrt(((normalized[:, None, :] - normalized[None, :, :]) ** 2).sum(axis=2))
    crowding = (distance < 0.13).sum(axis=1)
    order = sorted(range(len(data)), key=lambda i: (crowding[i], len(str(data.iloc[i]["city"]))), reverse=True)
    point_boxes = [plt.matplotlib.transforms.Bbox.from_bounds(px - 4, py - 4, 8, 8) for px, py in points]
    placed = []

    for i in order:
        row = data.iloc[i]
        label = str(row["city"])
        best = None
        for dx, dy in candidates:
            ha = "left" if dx > 1 else "right" if dx < -1 else "center"
            va = "bottom" if dy > 1 else "top" if dy < -1 else "center"
            trial = ax.annotate(
                label,
                xy=(row[x_col], row[y_col]),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=7.0,
                ha=ha,
                va=va,
            )
            box = trial.get_window_extent(renderer).expanded(1.05, 1.12)
            trial.remove()
            outside = box.width * box.height - _overlap_area(box, axes_box)
            label_overlap = sum(_overlap_area(box, other) for other in placed)
            point_overlap = sum(_overlap_area(box, other) for j, other in enumerate(point_boxes) if j != i)
            score = 1000.0 * (outside + label_overlap) + 20.0 * point_overlap + dx * dx + dy * dy
            if best is None or score < best[0]:
                best = (score, dx, dy, ha, va, box)

        _, dx, dy, ha, va, box = best
        ax.annotate(
            label,
            xy=(row[x_col], row[y_col]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.0,
            color="#222222",
            ha=ha,
            va=va,
            zorder=4,
            arrowprops={"arrowstyle": "-", "color": "#888888", "lw": 0.4, "alpha": 0.58, "shrinkA": 1.5, "shrinkB": 3.5},
        )
        placed.append(box)


def regression_ci(x: np.ndarray, y: np.ndarray, xx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    fit = stats.linregress(x, y)
    yy = fit.intercept + fit.slope * xx
    residual = y - (fit.intercept + fit.slope * x)
    dof = len(x) - 2
    s_err = np.sqrt(np.sum(residual**2) / dof)
    xbar = np.mean(x)
    sxx = np.sum((x - xbar) ** 2)
    se = s_err * np.sqrt((1 / len(x)) + ((xx - xbar) ** 2 / sxx))
    crit = stats.t.ppf(0.975, dof)
    return yy, yy - crit * se, yy + crit * se, fit.slope, fit.pvalue, fit.rvalue**2


def draw_quadrant(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    title: str,
    stem: str,
    night_greater_region: str,
) -> dict[str, float | int | str]:
    x = df[x_col].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)
    x_med = float(np.median(x))
    y_med = float(np.median(y))

    fig, ax = plt.subplots(figsize=(8.15, 6.55))
    axis_min = max(0.0, np.floor(min(float(x.min()), float(y.min())) * 10.0) / 10.0)
    axis_max = 1.00

    xx = np.linspace(axis_min, axis_max, 220)
    if night_greater_region == "above":
        ax.fill_between(xx, xx, axis_max, color="#EDEDED", alpha=0.75, linewidth=0, zorder=-2)
    else:
        ax.fill_between(xx, axis_min, xx, color="#EDEDED", alpha=0.75, linewidth=0, zorder=-2)

    yy, lo, hi, beta, p_value, r2 = regression_ci(x, y, xx)
    ax.fill_between(xx, lo, hi, color="#F2A37D", alpha=0.16, linewidth=0, zorder=-1)
    ax.plot(xx, yy, color="#D94B3D", lw=1.55, zorder=1.5)
    ax.plot([axis_min, axis_max], [axis_min, axis_max], color="#EFA0A0", lw=1.15, ls=(0, (1.5, 2.2)), zorder=0)
    ax.axvline(x_med, color="#AFAFAF", lw=1.0, ls=(0, (4, 3)), zorder=0)
    ax.axhline(y_med, color="#AFAFAF", lw=1.0, ls=(0, (4, 3)), zorder=0)

    for region, marker in MARKERS.items():
        reg = df[df["region"].eq(region)]
        for group, color in COLORS.items():
            sub = reg[reg["gii_group"].eq(group)]
            if sub.empty:
                continue
            ax.scatter(
                sub[x_col],
                sub[y_col],
                s=72,
                marker=marker,
                color=color,
                edgecolor="white",
                linewidth=0.75,
                alpha=0.92,
                zorder=3,
            )

    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    ax.set_aspect("equal", adjustable="box")
    ticks = np.arange(axis_min, 1.01, 0.1 if axis_min < 0.2 else 0.2)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel(x_label, fontsize=12.5, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=12.5, fontweight="bold")
    ax.set_title(title, loc="left", fontsize=14.5, fontweight="bold", pad=8)
    ax.tick_params(axis="both", labelsize=10, length=4, width=0.9)

    color_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="white", markeredgewidth=0.75, markersize=7.8, label=group.replace(" GII", ""))
        for group, color in COLORS.items()
    ]
    marker_handles = [
        Line2D([0], [0], marker=marker, color="#888888", linestyle="none", markerfacecolor="#BDBDBD", markeredgecolor="white", markeredgewidth=0.75, markersize=7.8, label=region)
        for region, marker in MARKERS.items()
    ]
    reference_handles = [
        Line2D([0], [0], color="#D94B3D", lw=1.55, label=f"Linear fit (R2 = {r2:.2f})"),
        Line2D([0], [0], color="#F2A37D", lw=6.0, alpha=0.18, label="95% confidence band"),
        Line2D([0], [0], color="#EFA0A0", lw=1.2, ls=(0, (1.5, 2.2)), label="Daytime = nighttime"),
        Line2D([0], [0], color="#D9D9D9", lw=6.0, alpha=0.75, label="Nighttime gap > daytime gap"),
        Line2D([0], [0], color="#AFAFAF", lw=1.0, ls=(0, (4, 3)), label="Sample medians"),
    ]

    leg1 = ax.legend(handles=color_handles, title="GII group", loc="upper left", bbox_to_anchor=(1.02, 0.97), fontsize=8.4, title_fontsize=9.0)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=marker_handles + reference_handles, title="Region / reference", loc="upper left", bbox_to_anchor=(1.02, 0.65), fontsize=8.1, title_fontsize=9.0)
    ax.add_artist(leg2)

    fig.tight_layout()
    annotate_cities(ax, df.reset_index(drop=True), x_col, y_col)
    for ext in ["png", "svg"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=450 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {
        "figure": stem,
        "x_col": x_col,
        "y_col": y_col,
        "n_cities": len(df),
        "x_median": x_med,
        "y_median": y_med,
        "beta": beta,
        "p_value": p_value,
        "r2": r2,
        "n_night_above_day": int((df["night_inequality"] > df["day_inequality"]).sum()),
    }


def main() -> None:
    apply_style()
    df = pd.read_csv(DATA)

    summaries = [
        draw_quadrant(
            df,
            x_col="day_inequality",
            y_col="night_inequality",
            x_label="Daytime gender gap",
            y_label="Nighttime gender gap",
            title="Gender inequality quadrant analysis",
            stem="gender_inequality_quadrant_gii_with_labels_day_x_night_y",
            night_greater_region="above",
        ),
        draw_quadrant(
            df,
            x_col="night_inequality",
            y_col="day_inequality",
            x_label="Nighttime gender gap",
            y_label="Daytime gender gap",
            title="Gender inequality quadrant analysis",
            stem="gender_inequality_quadrant_gii_with_labels_night_x_day_y",
            night_greater_region="below",
        ),
    ]

    summary = pd.DataFrame(summaries)
    summary.to_csv(LOG_DIR / "S1_gender_inequality_quadrant_summary.csv", index=False, encoding="utf-8-sig")
    manifest = "\n".join(p.relative_to(ROOT).as_posix() for p in sorted(FIG_DIR.glob("*")))
    (LOG_DIR / "figure_manifest.txt").write_text(manifest, encoding="utf-8")
    print(manifest)


if __name__ == "__main__":
    main()
