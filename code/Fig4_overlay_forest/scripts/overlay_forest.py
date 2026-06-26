from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
LOG_DIR = ROOT / "logs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

COL_FEMALE = "#7B2D8E"
COL_MALE = "#E6863B"
TEXT = "#222222"
GRID = "#BDBDBD"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["font.size"] = 10


ORDER = [
    ("NTL & Spatial\nStructure", "ntl_mean_log_z", "* Nightlight Intensity (NTL)"),
    ("NTL & Spatial\nStructure", "ntl_alpha_z", "NTL Spatial Decay (alpha)"),
    ("NTL & Spatial\nStructure", "ntl_D_z", "NTL Spatial Extent (D)"),
    ("NTL & Spatial\nStructure", "ntl_c_z", "NTL Core Brightness (c)"),
    ("NTL & Spatial\nStructure", "ulpi_z", "Urban Light-Population Index"),
    ("Built\nEnvironment", "walkability_global_z", "Walkability Index"),
    ("Built\nEnvironment", "ssx_choice_r5000_z", "Street Connectivity (SSX)"),
    ("Built\nEnvironment", "NDVI_mean_z", "Vegetation Cover (NDVI)"),
    ("Built\nEnvironment", "built_density_2020_z", "Building Density"),
    ("Population and\nDemographic", "Pop_Female_20_50_z", "Female Population (20-50 yrs)"),
    ("Population and\nDemographic", "ppl35prop_z", "Post-Childbearing Age (35+ Prop.)"),
    ("Socioeconomic\nconditions", "gdp_2022_z", "GDP per Grid Cell"),
    ("Gender related\nurban well-being", "Safety_Walk_Night_z", "Safety: Walking at Night"),
    ("Gender related\nurban well-being", "Female_Employment_Rate_15plus_z", "Female Employment Rate (15+)"),
    ("Gender related\nurban well-being", "SDG_3_1_Maternal_mortality_ratio_z", "SDG 3.1 Maternal Mortality"),
    ("Gender related\nurban well-being", "SDG_3_7_Adolescent_birth_rate_z", "SDG 3.7 Adolescent Birth Rate"),
    ("Gender related\nurban well-being", "SDG_4_4_Education_F_z", "SDG 4.4 Female Education"),
    ("Gender related\nurban well-being", "SDG_5_5_Share_of_seats_in_parliament_z", "SDG 5.5 Women in Parliament"),
    ("Gender related\nurban well-being", "SDG11_71pct_z", "SDG 11.7.1 Public Open Space"),
    ("Climate &\ntemporal effects", "tavg_night_z", "Night Temperature"),
    ("Climate &\ntemporal effects", "prcp_night_z", "Night Precipitation"),
    ("Climate &\ntemporal effects", "wspd_night_z", "Night Wind Speed"),
    ("Climate &\ntemporal effects", "pres_night_z", "Night Pressure"),
    ("Climate &\ntemporal effects", "rhum_night_z", "Night Humidity"),
    ("Climate &\ntemporal effects", "is_weekend", "Weekend Effect"),
]


def stars(p_value: float | None) -> str:
    if p_value is None or pd.isna(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def fmt_ci(row: pd.Series) -> str:
    sig = stars(row.get("p_value"))
    prefix = f"{sig} " if sig else ""
    return f"{prefix}{row['estimate']:+.3f} [{row['ci_low']:+.3f}, {row['ci_high']:+.3f}]"


def ordered_rows(df: pd.DataFrame, include_ntl: bool) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    terms = set(df["term"])
    for group, term, label in ORDER:
        if term == "ntl_mean_log_z" and not include_ntl:
            continue
        if term in terms:
            rows.append((group, term, label))
    return rows


def group_ranges(rows: list[tuple[str, str, str]]) -> list[tuple[str, float, float, float]]:
    ranges: list[tuple[str, int, int]] = []
    start = 0
    current = rows[0][0]
    for i, (group, _, _) in enumerate(rows):
        if group != current:
            ranges.append((current, start, i - 1))
            start = i
            current = group
    ranges.append((current, start, len(rows) - 1))
    out: list[tuple[str, float, float, float]] = []
    n = len(rows)
    for group, first, last in ranges:
        top = n - 1 - first
        bottom = n - 1 - last
        out.append((group, bottom, top, (top + bottom) / 2.0))
    return out


def draw_forest(df: pd.DataFrame, include_ntl: bool, stem: str, title: str, wide_x: bool = True) -> None:
    rows = ordered_rows(df, include_ntl=include_ntl)
    n = len(rows)
    y_map = {term: n - 1 - i for i, (_, term, _) in enumerate(rows)}
    df = df[df["term"].isin(y_map)].copy()
    df["y"] = df["term"].map(y_map)

    x_min = float(df["ci_low"].min())
    x_max = float(df["ci_high"].max())
    max_abs = max(abs(x_min), abs(x_max))
    lim = max(0.36, max_abs * (1.35 if wide_x else 1.12))
    width = 23.5 if wide_x else 16.0
    ratios = [1.55, 3.25, 9.9, 4.25] if wide_x else [1.55, 3.05, 4.8, 3.25]
    height = max(10.5, 0.42 * n + 2.0)

    fig = plt.figure(figsize=(width, height))
    gs = GridSpec(1, 4, width_ratios=ratios, wspace=0.035, figure=fig)
    ax_group = fig.add_subplot(gs[0, 0])
    ax_label = fig.add_subplot(gs[0, 1], sharey=ax_group)
    ax = fig.add_subplot(gs[0, 2], sharey=ax_group)
    ax_text = fig.add_subplot(gs[0, 3], sharey=ax_group)

    for axis in (ax_group, ax_label, ax_text):
        axis.set_ylim(-0.75, n - 0.25)
        axis.set_xlim(0, 1)
        axis.axis("off")

    ax.set_ylim(-0.75, n - 0.25)
    ax.set_xlim(-lim, lim)
    ax.axvline(0, color="#6F6F6F", linewidth=0.8, linestyle=(0, (5, 4)), zorder=0)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=9, colors=TEXT)
    ax.set_xlabel("Standardised coefficient (95% CI)", fontsize=10.5, fontweight="bold", labelpad=8)

    for y, (group, term, label) in zip(np.arange(n)[::-1], rows):
        ax_label.text(
            0.02,
            y,
            label,
            ha="left",
            va="center",
            fontsize=10.4,
            fontstyle="italic" if term == "ntl_mean_log_z" and include_ntl else "normal",
            fontweight="bold" if term == "ntl_mean_log_z" and include_ntl else "normal",
            color=TEXT,
        )

    for group, bottom, top, mid in group_ranges(rows):
        ax_group.text(
            0.94,
            mid,
            group,
            ha="right",
            va="center",
            fontsize=10.1,
            fontweight="bold",
            fontstyle="italic",
            linespacing=1.08,
            color="#111111",
        )
    for _, bottom, _, _ in group_ranges(rows)[1:]:
        sep_y = bottom - 0.5
        for axis in (ax_group, ax_label, ax, ax_text):
            axis.axhline(sep_y, color=GRID, linewidth=0.75, linestyle=(0, (4, 4)), zorder=0)

    offsets = {"Female PA": 0.13, "Male PA": -0.13}
    colors = {"Female PA": COL_FEMALE, "Male PA": COL_MALE}
    for outcome in ["Female PA", "Male PA"]:
        sub = df[df["outcome"] == outcome].copy()
        y = sub["y"].to_numpy(dtype=float) + offsets[outcome]
        x = sub["estimate"].to_numpy(dtype=float)
        xerr = np.vstack([x - sub["ci_low"].to_numpy(dtype=float), sub["ci_high"].to_numpy(dtype=float) - x])
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            fmt="o",
            markersize=4.0,
            linewidth=1.9 if wide_x else 1.35,
            capsize=0,
            color=colors[outcome],
            ecolor=colors[outcome],
            label=outcome,
            zorder=3,
        )

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.075), ncol=2, frameon=False, fontsize=10)
    ax_text.text(0.02, n - 0.02, "Female est. [95% CI]", ha="left", va="bottom", fontsize=9.5, fontweight="bold", color=COL_FEMALE)
    ax_text.text(0.52, n - 0.02, "Male est. [95% CI]", ha="left", va="bottom", fontsize=9.5, fontweight="bold", color=COL_MALE)
    for y, (_, term, _) in zip(np.arange(n)[::-1], rows):
        for outcome, x_pos in [("Female PA", 0.02), ("Male PA", 0.52)]:
            match = df[(df["term"] == term) & (df["outcome"] == outcome)]
            if match.empty:
                continue
            ax_text.text(x_pos, y, fmt_ci(match.iloc[0]), ha="left", va="center", fontsize=8.45, color=TEXT)

    fig.suptitle(title, x=0.02, y=0.992, ha="left", va="top", fontsize=13.5, fontweight="bold")
    fig.subplots_adjust(left=0.035, right=0.985, top=0.935, bottom=0.075, wspace=0.03)
    for ext in ["png", "svg"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=450 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def main() -> None:
    lmm = pd.read_csv(DATA_DIR / "lmm_plot_data.csv")
    gamm = pd.read_csv(DATA_DIR / "gamm_plot_data.csv")
    draw_forest(
        lmm,
        include_ntl=True,
        stem="lmm_complete_PA_overlay_right_label_columns_wide_x",
        title="LMM fixed-effect coefficients: female and male PA volume",
    )
    draw_forest(
        gamm,
        include_ntl=False,
        stem="gamm_full_PA_overlay_right_label_columns_wide_x",
        title="GAMM parametric coefficients: female and male PA volume",
    )
    manifest = "\n".join(p.relative_to(ROOT).as_posix() for p in sorted(FIG_DIR.glob("*")))
    (LOG_DIR / "figure_manifest.txt").write_text(manifest, encoding="utf-8")
    print(manifest)


if __name__ == "__main__":
    main()
