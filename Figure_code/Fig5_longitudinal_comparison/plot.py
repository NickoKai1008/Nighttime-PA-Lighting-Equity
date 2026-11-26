"""Draw Figure 5 from its panel-level source data."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Patch, Polygon, Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
WIDTH, HEIGHT = 1007, 1148  # Layout coordinates measured from the upper left
PURPLE, DARK_PURPLE, ORANGE, DARK_ORANGE = "#933b9e", "#732581", "#fb9b51", "#d46a1a"
TEAL, GREY, PALE, INK = "#147d79", "#85919c", "#dfe2e5", "#1c2a33"
GRID = "#eceef0"
DAY_NIGHT = LinearSegmentedColormap.from_list(
    "day_night", ["#fae9af", "#e5b16e", "#948894", "#474963", "#1c2038"])
font = next(name for name in ["Arial", "Liberation Sans", "DejaVu Sans"]
            if name in {f.name for f in font_manager.fontManager.ttflist})
plt.rcParams.update({
    "font.family": font, "font.size": 6.3, "axes.labelcolor": INK,
    "axes.edgecolor": "#53606e", "xtick.color": "#53606e", "ytick.color": "#53606e",
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "xtick.major.pad": 2, "ytick.major.pad": 2,
    "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42})


def table(name):
    return pd.read_csv(DATA / f"{name}.csv")


def axes(fig, x, y, width, height, **kwargs):
    return fig.add_axes([x / WIDTH, 1 - (y + height) / HEIGHT,
                         width / WIDTH, height / HEIGHT], **kwargs)


def label(fig, x, y, text, **kwargs):
    return fig.text(x / WIDTH, 1 - y / HEIGHT, text, va="center", **kwargs)


def grid(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.5)


def number(value, places=2):
    return f"{value:+.{places}f}".replace("-", "−")


def interval(low, high):
    separator = " to " if low < 0 else "–"
    return f"{low:.2f}{separator}{high:.2f}".replace("-", "−")


def panel_darkness(fig):
    days = table("a_dark_minutes")
    curves = table("a_darkness_curves")
    effects = table("a_darkness_estimates").set_index("unit")
    spline = table("a_radiance_curve")
    steps = table("a_radiance_steps").set_index("contrast")
    headings = [(94, "Natural darkness changes"),
                (346, "More dark minutes, stronger activity–radiance association"),
                (768, "…larger contrast on dimmer streets")]
    for x, title in headings:
        label(fig, x, 25, title, fontsize=6.15, weight="bold", color=INK)
    for start, end in [(270, 330), (696, 754)]:
        fig.add_artist(FancyArrowPatch((start / WIDTH, 1 - 27 / HEIGHT),
                                       (end / WIDTH, 1 - 27 / HEIGHT),
                                       transform=fig.transFigure, arrowstyle="-|>",
                                       mutation_scale=6, lw=0.55, color=PURPLE))
    for x, y, s in [(86, 44, "same streets, same clock hour, different days"),
                    (383, 47, "fixed radiance per cell; women's recorded PA"),
                    (720, 49, "same radiance step (25% of median)")]:
        label(fig, x, y, s, color=GREY, fontsize=5.05)
    ax = axes(fig, 80, 54, 239, 146)
    dark = days.dark_minutes_18_19
    ax.bar(days.doy, dark, width=1, color=DAY_NIGHT(dark / 60), lw=0)
    midmonth = pd.to_datetime([f"2025-{m:02d}-15" for m in range(1, 13)])
    ax.set(xticks=midmonth.dayofyear, xticklabels=list("JFMAMJJASOND"),
           xlim=(1, 366), ylim=(0, 60), yticks=[0, 30, 60],
           ylabel="Dark minutes\nevery 1 h in evening 17:00–\n22:00; 18:00–19:00 shown")
    ax.set_ylabel(ax.get_ylabel(), fontsize=5.9, labelpad=1.5)
    ax.tick_params(axis="x", length=0)
    grid(ax)
    ax.text(240, 52, "Melbourne 2025", fontsize=5.1)
    ax.text(4, 9, "Summer:\nhour in daylight", fontsize=5.25, style="italic", color=INK)
    ax.text(119, 35, "Winter:\nhour in darkness", fontsize=5.25, style="italic", color="white")
    icon = axes(fig, 96, 110, 35, 35)
    icon.set(xlim=(0, 1), ylim=(0, 1), aspect="equal")
    icon.axis("off")
    icon.add_patch(Circle((0.5, 0.5), 0.22, color="#fbb856"))
    for angle in np.arange(8) * np.pi / 4:
        r = np.array([0.33, 0.49])
        icon.plot(0.5 + r * np.cos(angle), 0.5 + r * np.sin(angle), color="#f2ac3e", lw=.8)
    icon = axes(fig, 178, 58, 30, 30)
    icon.set(xlim=(0, 1), ylim=(0, 1), aspect="equal")
    icon.axis("off")
    icon.add_patch(Circle((.48, .5), .31, color="#f5edd9"))
    icon.add_patch(Circle((.62, .61), .25, color="#1c2038"))
    icon.plot([.1, .92], [.88, .2], ".", color="#f5edd9", ms=1)

    ax = axes(fig, 370, 54, 287, 116)
    mel = curves[curves.unit == "melbourne_2025"]
    pool = curves[curves.unit == "core"]
    ax.fill_between(mel.dark_minutes, mel.lo, mel.hi, color="#e8d6ec")
    ax.plot(pool.dark_minutes, pool.pct, color=GREY, lw=1, ls=(0, (2, 1.5)))
    ax.plot(mel.dark_minutes, mel.pct, color=DARK_PURPLE, lw=1.65)
    ax.set(xlim=(0, 60), ylim=(0, 30), yticks=[0, 10, 20, 30], xticks=[],
           ylabel="Change in activity ratio per\ndoubling of radiance (%)")
    ax.set_ylabel(ax.get_ylabel(), fontsize=5.8, labelpad=1)
    grid(ax)
    row = effects.loc["melbourne_2025"]
    ax.text(1.5, 25.1, f"{number(row.estimate_percent)}%", color=DARK_PURPLE, fontsize=7.4, weight="bold")
    ax.text(1.5, 22.4, f"per 10 dark minutes\nWomen 95% CI {row.lower95_percent:.2f}–{row.upper95_percent:.2f}%",
            fontsize=4.8, color=INK, va="top")
    row = effects.loc["core"]
    ax.text(28, 14.2, f"Eight-city pool {number(row.estimate_percent)}% ({interval(row.lower95_percent, row.upper95_percent)})",
            fontsize=5.2, rotation=23, color=GREY)
    ax.text(44, 13.1, "Melbourne", color=DARK_PURPLE, fontsize=6.8, style="italic", weight="bold", rotation=11)
    bar = axes(fig, 370, 170, 287, 30)
    bar.imshow(np.linspace(0, 1, 256)[None, :], cmap=DAY_NIGHT, extent=(0, 60, 0, 1), aspect="auto")
    bar.set(yticks=[], xticks=[0, 15, 30, 45, 60], xlabel="Minutes of darkness within a fixed clock hour")
    bar.set_xlabel(bar.get_xlabel(), fontsize=5.4, labelpad=1)
    bar.tick_params(axis="x", labelsize=5.4, length=2)
    bar.text(1.5, .5, "More daylight minutes", va="center", fontsize=5.5, style="italic")
    bar.text(58.5, .5, "More dark minutes", va="center", ha="right", fontsize=5.5, style="italic", color="white")

    ax = axes(fig, 710, 55, 271, 145)
    subset = spline[spline.relative_radiance <= 6.12]
    ax.fill_between(subset.relative_radiance, subset.lower95_pct, subset.upper95_pct, color="#e8d6ec")
    ax.plot(subset.relative_radiance, subset.estimate_pct, color=DARK_PURPLE, lw=1.4)
    ax.axhline(0, color=GREY, lw=.5)
    ax.axvline(1, color=GREY, lw=.5, ls=":")
    ax.axhspan(-15, -7, color="#757887", zorder=0)
    for x, key, subtitle, xpos in [(.4, "F_low", "P25 - dimmer street", .4),
                                   (4.28, "F_high", "P75 - brighter street", 3.85)]:
        step = steps.loc[key]
        low, high = step.start_relative_radiance, step.end_relative_radiance
        ax.axvspan(low, high, color=TEAL, alpha=.30, lw=0)
        radiance = np.r_[low, subset.relative_radiance[subset.relative_radiance.between(low, high)], high]
        fitted = np.interp(np.log(radiance), np.log(spline.relative_radiance), spline.estimate_pct)
        ax.plot(radiance, fitted, color=DARK_PURPLE, lw=2.7, solid_capstyle="round")
        ax.text(x, 8.1, f"{number(step.pct)}%", color=DARK_PURPLE, weight="bold", fontsize=7.3)
        ax.text(xpos, 5.5, subtitle, color=INK, fontsize=5.6)
    ax.text(1.06, -5.6, "city median", color=GREY, fontsize=5.4)
    for x, opacity in [(.43, .3), (3.25, .7), (3.85, .7), (4.45, .7)]:
        ax.add_patch(Polygon([(x, -9.15), (x - .32, -15), (x + .32, -15)],
                             closed=True, fc="#e6b54f", ec="none", alpha=opacity))
        ax.plot([x, x], [-15, -9.15], color="#ccd1d1", lw=.6)
    ax.text(.49, -10.4, "Dim", fontsize=5.8, style="italic", color="white")
    ax.text(4.69, -10.4, "Bright", fontsize=5.8, style="italic", color="white")
    ax.set(xlim=(0, 6), ylim=(-15, 12), yticks=[-10, 0, 10], xticks=[0, 1, 2, 4, 6],
           ylabel="Darkness association (%)", xlabel="Street radiance relative to city median")
    ax.set_xlabel(ax.get_xlabel(), labelpad=1, fontsize=5.5)
    ax.tick_params(labelsize=5.7)


def panel_radiance(fig):
    quarterly = table("b_radiance_quarterly")
    dates = pd.to_datetime(quarterly.quarter)
    ax = axes(fig, 80, 268, 582, 232)
    ax.fill_between(dates, quarterly.ci_low, quarterly.ci_high, lw=0, color="#e7e9eb")
    ax.plot(dates, quarterly.estimate, color="#687279", lw=1.15)
    for row in table("b_radiance_annual").itertuples():
        date = pd.Timestamp(row.year, 7, 1)
        col = TEAL if row.year >= 2023 else (INK if row.year <= 2019 else "#96a3ac")
        ax.plot([date, date], [row.ci_low, row.ci_high], color=col, lw=.75)
        ax.plot(date, row.estimate, "o", ms=4.1, mfc=col, mec="white", mew=.6)
    windows = {"pre": pd.Timestamp("2019-11-01"), "post": pd.Timestamp("2023-11-01")}
    for text, date in windows.items():
        align = "left" if text == "pre" else "right"
        ax.axvline(date, color=TEAL, lw=.8, ls="--")
        ax.annotate(text, (date, 29), xytext=((2 if align == "left" else -4), 0),
                    textcoords="offset points", ha=align, color=TEAL)
    row = table("b_radiance_contrast").iloc[0]
    ax.text(pd.Timestamp("2021-11-01"), 27, f"DID {number(row.estimate)}\n({interval(row.ci_low, row.ci_high)})",
            ha="center", va="top", color=INK, weight="bold", fontsize=7.0)
    ax.axhline(0, color=GREY, lw=.5)
    years = np.arange(2017, 2027)
    ax.set(ylim=(-17, 31), yticks=np.arange(-16, 25, 8), xlim=(pd.Timestamp("2016-10-01"), pd.Timestamp("2026-01-15")),
           xticks=pd.to_datetime([f"{y}-01-01" for y in years]), xticklabels=[str(y) if y % 2 else "" for y in years],
           ylabel="Treatment − control radiance\n(nW cm$^{-2}$ sr$^{-1}$)")
    grid(ax)
    boxes = table("b_image_annotations")
    pair = table("b_image_pair").iloc[0]
    for x, period, filename in [(178, "pre", pair.before_image), (529, "post", pair.after_image)]:
        photo = axes(fig, x, 238, 88, 93)
        with Image.open(DATA / filename) as im:
            photo.imshow(np.asarray(im), extent=(0, 800, 800, 0), aspect="auto")
        for box in boxes[boxes.period == period].itertuples():
            context = box.role == "context"
            photo.add_patch(Rectangle((box.x1, box.y1), box.x2 - box.x1, box.y2 - box.y1,
                                      fill=False, ec=GREY if context else "#f4cf55", lw=.45 if context else .75))
        photo.set(xlim=(0, 800), ylim=(800, 0), xticks=[], yticks=[])
        for spine in photo.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(.55)
        photo.text(.03, .05, f"{windows[period].year} {period}", transform=photo.transAxes,
                   weight="bold", color="white", fontsize=6.5,
                   bbox={"fc": INK, "ec": "none", "alpha": .55, "pad": .8})
    ring = axes(fig, 551, 419, 79, 80)
    ring.set(xlim=(-1.6, 1.6), ylim=(-1.6, 1.6), aspect="equal")
    ring.axis("off")
    for radius, col in [(1.25, PALE), (.83, "white"), (.21, TEAL)]:
        ring.add_patch(Circle((0, 0), radius, color=col, lw=0))
    ring.text(0, .82, "Control", ha="center", color=INK)
    ring.text(0, -.66, "Treatment", ha="center", color=TEAL)


def panel_streetscape(fig):
    labels = {"Road surface": "Road", "Walking surfaces": "Walkway", "Tree canopy": "Trees",
              "Low vegetation": "Low veg.", "Built frontage": "Buildings", "Sky": "Sky",
              "Water": "Water", "Barriers": "Barriers", "Lighting fixtures": "Lighting",
              "Poles": "Poles", "Street furniture": "Furniture", "Vehicles & people": "Vehicles,\npeople"}
    rows = table("c_streetscape").set_index("group").loc[list(labels)]
    angles = np.linspace(0, 2 * np.pi, len(rows), endpoint=False)
    ax = axes(fig, 713, 289, 225, 225, projection="polar")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    for col, color, name in [("before_pct", GREY, "Before"), ("after_pct", TEAL, "After")]:
        values = rows[col].to_numpy()
        theta, radius = np.r_[angles, angles[0]], np.r_[values, values[0]]
        ax.plot(theta, radius, color=color, lw=1.35, label=name)
        ax.fill(theta, radius, color=color, alpha=.17, lw=0)
    ax.set(ylim=(0, 40), yticks=[10, 20, 30, 40], yticklabels=["", "20%", "", "40%"],
           xticks=angles, xticklabels=list(labels.values()))
    ax.set_rlabel_position(165)
    ax.tick_params(axis="x", pad=-1, labelsize=6)
    ax.tick_params(axis="y", labelsize=5.7)
    ax.grid(color=PALE, lw=.5)
    ax.spines["polar"].set_color(PALE)
    ax.legend(loc="lower right", bbox_to_anchor=(1.22, -.10), fontsize=6.1, handlelength=1.2)


def panel_day_night(fig):
    projects, estimates = table("d_project_changes"), table("d_group_changes")
    ax = axes(fig, 82, 560, 262, 160)
    for start, outcome, light, dark, sex in [(0, "female_pa", "#b66ac4", DARK_PURPLE, "Women"),
                                              (2.25, "male_pa", "#ffac71", DARK_ORANGE, "Men")]:
        xs = [start, start + 1]
        pairs = projects[projects.measure == outcome].pivot(index="project_id", columns="period", values="did")
        for pair in pairs.itertuples():
            ax.plot(xs, [pair.day, pair.night], "-o", color=light, alpha=.75, lw=.5, ms=2, mew=0)
        values = estimates[estimates.outcome == outcome].set_index("estimand").loc[["day_did", "night_did"], "estimate"]
        ax.plot(xs, values, "-D", lw=2.3, color=dark, mec="white", mew=.8, ms=5.3)
        for i, y in enumerate(values):
            ax.annotate(number(y, 0), (xs[i], y), xytext=((-4 if i == 0 else 4), 0), textcoords="offset points",
                        ha="right" if i == 0 else "left", va="center", color=dark, weight="bold", fontsize=6.4)
        ax.text(start + .5, 1, sex, transform=ax.get_xaxis_transform(), va="bottom", ha="center",
                weight="bold", color=dark, fontsize=6.4)
    ax.set_yscale("symlog", linthresh=50, linscale=.9)
    ax.axhline(0, color=GREY, lw=.5)
    ax.set(ylim=(-700, 2600), xlim=(-.75, 4), xticks=[0, 1, 2.25, 3.25],
           xticklabels=["Day", "Night", "Day", "Night"], yticks=[-500, -50, 0, 50, 500, 2000],
           yticklabels=["−500", "−50", "0", "50", "500", "2,000"],
           ylabel="Project − reference change\n(records per 200-m cell)")
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    grid(ax)


def panel_placebos(fig, letter):
    estimates = table("ef_estimates").set_index(["outcome", "estimand"])
    draws = table("ef_spatial_placebos")
    if letter == "e":
        ax = axes(fig, 396, 564, 238, 109)
        groups = [("Women", "female_pa", "ddd", PURPLE), ("Men", "male_pa", "ddd", ORANGE)]
        limit, ticks = 200, [-200, -100, 0, 100, 200]
        xlabel = "Night − day DDD (records per 200-m cell)"
    else:
        ax = axes(fig, 740, 564, 201, 109)
        groups = [("Night", "female_share", "night_did", PURPLE), ("Day", "female_share", "day_did", "#bf96c9")]
        limit, ticks = 20, [-20, -10, 0, 10, 20]
        xlabel = "DID (percentage points)"
    xs = np.linspace(-limit, limit, 600)
    for i, (name, outcome, estimand, color) in enumerate(groups):
        group = draws[(draws.outcome == outcome) & (draws.estimand == estimand)]
        density = gaussian_kde(group.placebo_estimate)(xs)
        ridge = i + .06 + .38 * density / density.max()
        ax.fill_between(xs, i + .06, ridge, lw=0, color=PALE)
        ax.plot(xs, ridge, color=GREY, lw=.5)
        row = estimates.loc[(outcome, estimand)]
        ax.plot([row.ci_low, row.ci_high], [i - .10, i - .10], lw=1.5, color=color)
        ax.plot(row.estimate, i - .1, "o", ms=4.7, color=color, mec="white", mew=.7)
        ax.annotate(name, (-limit, i - .1), xytext=(-2, 0), textcoords="offset points", ha="right", va="center")
        ax.annotate(number(row.estimate), (limit, i - .1), xytext=(3, 0), textcoords="offset points",
                    va="center", weight="bold" if i == 0 else "normal")
    ax.axvline(0, lw=.5, color=GREY, zorder=0)
    ax.spines["left"].set_visible(False)
    ax.set(xlim=(-limit, limit), ylim=(1.5, -.55), yticks=[], xticks=ticks, xlabel=xlabel)
    if letter == "e":
        ax.legend(handles=[Patch(fc=PALE, ec=GREY, lw=.4, label=f"Placebo areas (n = {len(group):,})")],
                  loc="upper right", bbox_to_anchor=(1.26, 1.15), handlelength=.65, fontsize=6.2)


CITIES = {"Hong Kong": ("o", "#1f9e89", 4.4), "Seoul": ("s", "#2f5fa7", 4.0),
          "Melbourne": ("^", "#d4aa2a", 4.8), "Shanghai": ("D", "#9a4f3a", 3.7)}


def scaled_change(percent):
    return np.log1p(np.asarray(percent) / 100) / np.log(1.25)


def panel_projects(fig):
    projects = table("g_projects")
    lines = table("g_fitted_lines")
    slopes = table("g_slopes").set_index("sex")
    diff = slopes.loc["F_minus_M"]
    label(fig, 154, 761, f"{int(diff.n_projects)} renewal projects in four cities", fontsize=6.6,
          color=INK, weight="bold", style="italic")
    label(fig, 366, 761, f"Women − men slope {number(diff.coefficient)} ({interval(diff.ci_low, diff.ci_high)})",
          color="#495763", fontsize=6.5, weight="bold", style="italic")
    left, right = scaled_change([-55, 48])
    for x, sex, name, color, light in [(60, "F", "Women", DARK_PURPLE, PURPLE), (364, "M", "Men", DARK_ORANGE, ORANGE)]:
        ax = axes(fig, x, 799, 285, 287)
        ax.add_patch(Rectangle((0, 0), right, 3.62, color=light, alpha=.07, lw=0))
        ax.add_patch(Rectangle((left, -3.62), -left, 3.62, color=light, alpha=.07, lw=0))
        rail = 4.05
        for y in (-rail, rail):
            ax.axhspan(y - .3, y + .3, color="#f1f2f3", lw=0)
        ax.axhline(0, color=GREY, lw=.6)
        ax.axvline(0, color=GREY, lw=.6)
        line = lines[lines.sex == sex]
        ax.fill_between(line.D_relative25, line.lo, line.hi, color=light, alpha=.20, lw=0)
        ax.plot(line.D_relative25, line.fit, color=color, lw=1.9)
        subset = projects[projects.sex == sex]
        for city, (marker, fill, size) in CITIES.items():
            rows = subset[subset.city_label == city]
            finite = rows[rows.ratio_status == "finite"]
            nonfinite = rows[rows.ratio_status != "finite"]
            ax.plot(finite.D_relative25, finite.adjusted_log_ddd.clip(-3.5, 3.5), marker, ms=size,
                    mfc=fill, mec="white", mew=.6, ls="none")
            ax.plot(nonfinite.D_relative25, np.sign(nonfinite.direction) * rail, marker, ms=size * .92,
                    mfc="white", mec=fill, mew=.85, ls="none")
        counts = subset.quadrant.value_counts()
        for quad, text, xx, yy, ha, emphasis, offset in [
            ("I", "Brighter, more at night", .98, .86, "right", True, -.075),
            ("II", "Dimmer, more at night", .02, .86, "left", False, -.075),
            ("III", "Dimmer, less at night", .02, .40, "left", True, .075),
            ("IV", "Brighter, less at night", .98, .13, "right", False, .075)]:
            ink = color if emphasis else GREY
            ax.text(xx, yy, text, transform=ax.transAxes, ha=ha, va="center", color=ink, fontsize=6.0,
                    weight="bold" if emphasis else "normal", bbox={"fc": "white", "ec": "none", "alpha": .75, "pad": .2})
            ax.text(xx, yy + offset, f"n = {counts.get(quad, 0)}", transform=ax.transAxes,
                    ha=ha, va="center", color=ink, fontsize=5.9)
        row = slopes.loc[sex]
        ax.text(0, 1.03, name, transform=ax.transAxes, color=color, weight="bold")
        ax.text(1, 1.03, f"slope {number(row.coefficient)} ({interval(row.ci_low, row.ci_high)})",
                transform=ax.transAxes, ha="right", color=INK)
        xt = [-40, -20, 0, 20, 40]
        ax.set(xlim=(left, right), ylim=(-4.55, 4.55), xticks=scaled_change(xt),
               xticklabels=[number(v, 0) if v else "0" for v in xt],
               yticks=[-rail, *np.log([1/16, 1/4, 1, 4, 16]), rail], xlabel="Radiance change vs reference (%)")
        if sex == "F":
            ax.set_yticklabels(["−∞", "1/16", "1/4", "1", "4", "16", "+∞"])
            ax.set_ylabel("Night-specific change, ratio to city expectation")
        else:
            ax.tick_params(labelleft=False)
    handles = [Line2D([], [], marker=mk, color="none", ms=sz + .4, mfc=col, mec="white", label=city)
               for city, (mk, col, sz) in CITIES.items()]
    handles.append(Line2D([], [], marker="o", ls="none", ms=4, mfc="white", mec=GREY, label="zero records in one cell"))
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.292, .011), ncol=5, fontsize=6.3,
               handlelength=1.1, columnspacing=1.25, handletextpad=.3)


def panel_starting_radiance(fig):
    curves = table("h_slope_curve")
    row = table("h_interaction").iloc[0]
    starts = table("h_starting_radiance").treatment_pre
    ax = axes(fig, 715, 799, 275, 287)
    ax.fill_between(curves.baseline_value, curves.ci_low, curves.ci_high, color=TEAL, alpha=.15, lw=0)
    ax.plot(curves.baseline_value, curves.coefficient, color=TEAL, lw=1.8)
    ax.axhline(0, color=GREY, lw=.65)
    rug = starts[starts.between(curves.baseline_value.min(), curves.baseline_value.max())]
    ax.plot(rug, np.full(len(rug), -1.6), "|", color=INK, ms=4, mew=.5)
    ax.text(.97, .77, f"{number(row.coefficient)} per doubling\n({interval(row.ci_low, row.ci_high)})",
            transform=ax.transAxes, ha="right", va="top", color=INK, linespacing=1.5)
    ax.text(.28, .56, "Smaller where already bright", transform=ax.transAxes,
            color=INK, weight="bold", style="italic", rotation=-17, fontsize=6.5)
    ax.set_xscale("log", base=2)
    ax.set(xticks=[50, 100, 150], xticklabels=["50", "100", "150"], ylim=(-1.85, 4.2), yticks=[-1, 0, 1, 2, 3, 4],
           xlabel="Starting radiance (nW cm$^{-2}$ sr$^{-1}$)", ylabel="Women's slope (per 25% radiance growth)")
    ax.minorticks_off()
    grid(ax)


def main():
    fig = plt.figure(figsize=(7.2, 7.2 * HEIGHT / WIDTH))
    panel_darkness(fig)
    panel_radiance(fig)
    panel_streetscape(fig)
    panel_day_night(fig)
    panel_placebos(fig, "e")
    panel_placebos(fig, "f")
    panel_projects(fig)
    panel_starting_radiance(fig)
    for letter, x, y in [("a", 22, 23), ("b", 22, 240), ("c", 702, 251), ("d", 22, 549),
                         ("e", 368, 549), ("f", 705, 549), ("g", 22, 763), ("h", 705, 763)]:
        label(fig, x, y, letter, fontsize=12.7, color=INK, weight="bold")
    out = ROOT / "output"
    out.mkdir(exist_ok=True)
    for extension in ["png", "svg", "pdf"]:
        fig.savefig(out / f"fig5.{extension}", dpi=600, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
