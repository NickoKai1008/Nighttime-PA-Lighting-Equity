"""Draw Figure 3 from the archived nighttime prediction and summary tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from matplotlib.ticker import NullLocator
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SEX = {"female": "#a442c6", "male": "#ff824b"}
LAND = {"recreational": "#43a568", "residential": "#377eb8", "commercial": "#fb792b"}
RAD = "NTL radiance (10⁻³ W m⁻² sr⁻¹ μm⁻¹)"
UNIT = "W m⁻² sr⁻¹ μm⁻¹"
XMAX = 17.453984458769434  # Pooled 99th percentile, in display units.
CITIES = ["hong_kong", "dubai", "sydney", "auckland", "kuala_lumpur",
          "doha", "bangkok", "tehran", "dhaka"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": .55, "xtick.major.width": .55,
                     "ytick.major.width": .55, "xtick.major.size": 2,
                     "ytick.major.size": 2, "xtick.major.pad": 2,
                     "ytick.major.pad": 2, "legend.frameon": False,
                     "svg.fonttype": "none", "pdf.fonttype": 42,
                     "svg.hashsalt": "ncities-figure-3", "savefig.facecolor": "white"})


def data(name):
    return pd.read_csv(ROOT / "data" / (name + ".csv"), float_precision="round_trip")


def axes(fig, x, y, width, height, **kwargs):
    """Place an axis in the manuscript layout, measured from its top left."""
    return fig.add_axes([x / 1133, 1 - (y + height) / 1500,
                         width / 1133, height / 1500], **kwargs)


def text_at(fig, x, y, label, **kwargs):
    return fig.text(x / 1133, 1 - y / 1500, label, va="center", **kwargs)


def ticks(axis, values):
    axis.set_ticks(values, labels=[f"{v:g}" for v in values])
    axis.set_minor_locator(NullLocator())


def pooled_panel(fig, curves, bands, bounds, hist, overlay):
    x, width = (563, 529) if overlay else (85, 392)
    ax = axes(fig, x, 35, width, 354)
    ax.set(xlim=(0, XMAX), yscale="log")
    scopes = list(LAND) if overlay else ["all25"]
    display_log = [0]
    for scope in scopes:
        for sex in ["male", "female"]:
            color = LAND[scope] if overlay else SEX[sex]
            p = curves[(curves.scope == scope) & (curves.sex == sex)].sort_values("L")
            b = bands[(bands.scope == scope) & (bands.sex == sex)].sort_values("L")
            for frame, fields in [(b, ["lo", "hi"]), (p, ["f"])]:
                radiance = frame.L.to_numpy() * 1000
                end = min(XMAX, radiance[-1])
                grid = np.r_[radiance[radiance < end], end]
                values = [np.interp(grid, radiance, frame[f]) for f in fields]
                display_log.extend([v.min() for v in values] + [v.max() for v in values])
                if len(values) == 2:
                    ax.fill_between(grid, *np.exp(values), color=color, edgecolor="none",
                                    alpha=.11 if sex == "male" else .17, zorder=1)
                else:
                    ax.plot(grid, np.exp(values[0]), color=color, lw=1,
                            ls="--" if overlay and sex == "male" else "-", zorder=3)
    ax.set_ylim(np.exp(min(display_log) - .15), np.exp(max(display_log) + .17))
    ticks(ax.yaxis, [v for v in [.25, .5, 1, 2, 4, 8, 16] if ax.get_ylim()[0] < v < ax.get_ylim()[1]])
    ax.axhline(1, color="#aaa", lw=.5, zorder=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for i, scope in enumerate(scopes):
        color = LAND[scope] if overlay else SEX["female"]
        row = bounds.loc[scope]
        xf = ax.get_xaxis_transform()
        for val, style, marker in [(row.onset, ":", "o"), (row.L90, "-", "|")]:
            at = val * 1000
            y = .07 + i * .075
            if at <= XMAX:
                ax.axvline(at, color=color, lw=.6, ls=style, zorder=2)
                ax.plot(at, y, marker=marker, ms=2.8, color=color, transform=xf)
            else:
                ax.plot(XMAX, y, marker=">", ms=3.5, color=color, mfc="white",
                        transform=xf, clip_on=False)
        if overlay:
            ax.text(row.onset * 1000 + .22, .071 + i * .095,
                    f"{row.onset:.4f}\n({UNIT})", transform=xf, fontsize=6.3,
                    va="center", zorder=4)
    hx = axes(fig, x, 406, width, 39, sharex=ax)
    hh = hist[hist.scope.isin(scopes)].pivot(index=["left", "right"], columns="scope",
                                             values="count_grid_days")[scopes]
    left = 1000 * hh.index.get_level_values("left").to_numpy()
    widths = 1000 * hh.index.get_level_values("right").to_numpy() - left
    bottom = np.zeros(len(hh))
    for scope in scopes:
        values = hh[scope].to_numpy() / hh.to_numpy().sum()
        hx.bar(left, values, bottom=bottom, width=widths, align="edge",
               color=LAND.get(scope, "#bbb"), alpha=.50, lw=0)
        bottom += values
    hx.spines["left"].set_visible(False)
    hx.set(yticks=[], xticks=[0, 4, 8, 12, 16], xlabel=RAD, xlim=(0, XMAX))
    hx.text(.985, .45, "Grid-day light distribution", ha="right",
            transform=hx.transAxes, fontsize=6.8, color="#555")
    return ax


def top_panels(fig):
    points = data("pooled_curves")
    bands = data("bootstrap_bands")
    bounds = data("female_thresholds").set_index("scope")
    hist = data("grid_day_histograms")
    ax = pooled_panel(fig, points, bands, bounds, hist, False)
    bx = pooled_panel(fig, points, bands, bounds, hist, True)
    ax.set_ylabel("Relative nighttime PA (×)", fontsize=8.5)
    limits = bounds.loc["all25", ["onset", "L90"]].to_numpy(float) * 1000
    ax.axvspan(*limits, color="#4aaaa4", alpha=.10, lw=0, zorder=0)
    ax.axvspan(limits[1], XMAX, color="#777777", alpha=.16, lw=0, zorder=0)
    ax.text(.04, .98, "Steeper association", va="top", ha="center", rotation=90,
            style="italic", weight="bold", transform=ax.transAxes)
    ax.text(.16, .98, "Diminishing\nassociation", va="top", ha="center", rotation=90,
            style="italic", weight="bold", transform=ax.transAxes)
    ax.text(.71, .98, "Flatter association", va="top", ha="center", rotation=90,
            style="italic", weight="bold", fontsize=6.8, transform=ax.transAxes)
    ax.text(.82, .38, "Energy saving\n&\nlight pollution\nmitigation potential",
            ha="center", va="center", transform=ax.transAxes, fontsize=6.8)
    female = points[(points.scope == "all25") & (points.sex == "female")].sort_values("L")
    photo_specs = [("low", 185, 227, limits[0]), ("mid", 296, 69, limits[1]),
                   ("high", 435, 175, 15.6)]
    for level, x, y, radiance in photo_specs:
        response = np.interp(radiance, female.L * 1000, np.exp(female.f))
        color = SEX["female"] if level == "mid" else "black"
        ax.plot(radiance, response, "o", ms=3.2, color=color, zorder=5)
        ax.annotate("", xy=(radiance, response), xytext=(x / 1133, 1 - y / 1500),
                    textcoords=fig.transFigure, arrowprops={"arrowstyle": "-", "lw": .55,
                    "ls": "--", "color": color})
        photo = axes(fig, x - 43, y - 43, 86, 86)
        im = photo.imshow(plt.imread(ROOT / "images" / f"street_{level}.jpg"))
        circle = Circle((.5, .5), .493, transform=photo.transAxes, fc="none", ec=color, lw=.6)
        im.set_clip_path(circle)
        photo.add_patch(circle)
        photo.set_axis_off()
    for radiance, tx, ty in [(limits[0], .12, .63), (limits[1], .53, .59)]:
        response = np.interp(radiance, female.L * 1000, np.exp(female.f))
        ax.annotate(f"{radiance / 1000:.4f}\n({UNIT})", xy=(radiance, response),
                    xytext=(tx, ty), textcoords="axes fraction", va="center",
                    ha="left" if tx < .5 else "center", fontsize=6.2,
                    arrowprops={"arrowstyle": "-", "lw": .5, "color": "#333"})
    inset = axes(fig, 148, 277, 312, 77, facecolor="none")
    ends = np.arange(2, 18, 2)
    for sex, offset in [("female", -.17), ("male", .17)]:
        g = points[(points.scope == "all25") & (points.sex == sex)].sort_values("L")
        ratios = np.interp(ends, g.L * 1000, np.exp(g.f))
        inset.bar(np.arange(7) + offset, np.diff(ratios), width=.26, color=SEX[sex])
    inset.set(xticks=np.arange(7), xticklabels=[f"{a}–{b}" for a, b in zip(ends[:-1], ends[1:])],
              yticks=[0, .2, .4], ylim=(0, .55))
    inset.tick_params(labelsize=5.1, length=1.5, pad=1)
    inset.set_ylabel("PA ratio change (ΔR)", fontsize=5.9, labelpad=0)
    inset.set_xlabel("Road-light interval (10⁻³ W m⁻² sr⁻¹ μm⁻¹)", fontsize=5.8, labelpad=1)
    inset.text(.10, .83, "Gain per equal light step", weight="bold", style="italic",
               fontsize=6.5, transform=inset.transAxes)
    handles = [Line2D([], [], color=color, lw=1.4) for color in SEX.values()]
    fig.legend(handles, ["Female", "Male"], loc="center left", bbox_to_anchor=(.10, .728),
               ncol=2, fontsize=6.8, handlelength=1.4, columnspacing=1.0)
    for scope, x in zip(LAND, [628, 792, 963]):
        text_at(fig, x, 38, "━━", color=LAND[scope])
        text_at(fig, x + 32, 38, scope.title(), fontsize=7.2)
        text_at(fig, x + 32, 62, f"{bounds.loc[scope, 'L90']:.4f}\n({UNIT})", fontsize=6.2)
    bx.legend([Line2D([], [], color="#333", ls=ls, lw=1) for ls in ["-", "--"]],
              ["Female", "Male"], loc="upper right", bbox_to_anchor=(1, .82), fontsize=7)
    bx.text(.56, .28, "│  L90 (90% response gain)", transform=bx.transAxes, fontsize=7)
    bx.text(.56, .15, "┊   Flattening onset", transform=bx.transAxes, fontsize=7)
    r = bands[(bands.scope == "recreational") & (bands.sex == "female")].iloc[0]
    bx.text(.99, .035, f"Accepted FRW: Recreational F {r.n_accepted}/{r.n_attempted}",
            ha="right", color="#555", transform=bx.transAxes, fontsize=6.5)


def context_panels(fig):
    scopes = ["all25", *LAND]
    ratios = data("onset_ratios").set_index("scope")
    ax = axes(fig, 160, 546, 267, 200)
    ax.axvline(1, color="#aaa", lw=.5)
    for y, scope in enumerate(reversed(scopes)):
        r = ratios.loc[scope]
        median = r.female_over_male_onset
        ax.errorbar(median, y, xerr=[[median - r.lo], [r.hi - median]],
                    fmt="o", color=SEX["female"], ms=3.3, lw=.8, capsize=1.6)
        ax.text(1.02, y, f"{median:.2f} [{r.lo:.2f}, {r.hi:.2f}]",
                transform=ax.get_yaxis_transform(), va="center", fontsize=7)
    ax.set(xscale="log", xlim=(.47, 9.8), ylim=(-.5, 3.5), yticks=range(4),
           yticklabels=["Commercial", "Residential", "Recreational", "All 25 cities"],
           xlabel="Female / male onset\nPaired FRW median (95% interval)")
    ticks(ax.xaxis, [.5, 1, 2, 4, 8])
    ax.set_title("Gender difference in onset", weight="bold", fontsize=8.2, pad=4)
    shares = data("exposure_shares").set_index(["scope", "weight"])
    bx = axes(fig, 765, 546, 334, 200)
    colors = ["#8e559e", "#40aaa6", "#ff8e56"]
    labels = []
    for i, scope in enumerate([*LAND, "three_land_uses"]):
        for j, weight in enumerate(["grid_days", "radiance"]):
            left = 0
            y = 7 - 2 * i - j
            row = shares.loc[scope, weight]
            for key, color in zip(["below", "within", "beyond"], colors):
                pct = 100 * row[key]
                bx.barh(y, pct, left=left, height=.72, color=color, lw=0, alpha=1 if j == 0 else .63)
                if pct > 8:
                    bx.text(left + pct / 2, y, f"{pct:.0f}%", ha="center", va="center",
                            color="white", weight="bold", fontsize=6.4)
                left += pct
            label = "All three" if i == 3 else scope.title()
            labels.append(label + (" · grid-days" if j == 0 else " · radiance-weighted"))
    bx.set(yticks=range(7, -1, -1), yticklabels=labels, ylim=(-.65, 7.65), xlim=(0, 100),
           xticks=[0, 25, 50, 75, 100], xticklabels=["0", "25", "50", "75", "100%"],
           xlabel="Share relative to female onset–L90")
    bx.tick_params(axis="y", length=0, labelsize=6.4)
    bx.set_title("Where grid-days and radiance fall", fontsize=8.2, weight="bold", pad=4)
    for x, label, color in zip([751, 866, 1007], ["Below onset", "Within range", "Beyond L90"], colors):
        text_at(fig, x, 798, label, color=color, fontsize=6.8)


def city_panels(fig):
    curves = data("city_curves")
    thresholds = data("city_thresholds").set_index("city")
    label_heights = [(.96, .78), (.68, .72), (.96, .82), (.96, .96),
                     (.96, .62), (), (.78,), (.38,), (.96, .76)]
    for index, city in enumerate(CITIES):
        x = [96, 439, 782][index % 3]
        y = [848, 1056, 1264][index // 3]
        ax = axes(fig, x, y, 286, 156)
        ax.set(xlim=(.08, 12), xscale="log", yscale="log")
        visible = [1.]
        for sex in ["male", "female"]:
            part = curves[(curves.city == city) & (curves.sex == sex)].sort_values("L")
            radiance = part.L.to_numpy() * 1000
            logged = part[["f", "lo", "hi"]].to_numpy()
            ratio, lower, upper = np.exp(logged).T
            visible.extend(np.exp(logged[radiance <= 12]).ravel())
            color = SEX[sex]
            ax.fill_between(radiance, lower, upper, color=color, lw=0,
                            alpha=.10 if sex == "male" else .17)
            ax.plot(radiance, ratio, lw=1.15, color=color)
            if radiance[-1] > 12:
                edge = [np.exp(np.interp(np.log(12), np.log(radiance), v)) for v in logged.T]
                visible.extend(edge)
                ax.plot(12, edge[0], ">", ms=2.7, color=color, clip_on=False)
        low, high = min(visible) * .84, max(visible) * 1.18
        ax.set_ylim(low, high)
        yt = np.array(sorted(v * 10.**power for v in [1, 2, 5] for power in range(-4, 4)))
        yt = yt[(yt >= low) & (yt <= high)]
        ticks(ax.yaxis, yt[::max(1, int(np.ceil(len(yt) / 5)))])
        ticks(ax.xaxis, [.1, 1, 10])
        ax.axhline(1, color="#afaab0", lw=.5, zorder=0)
        ax.set_title(city.replace("_", " ").title(), loc="left", fontsize=9.3, weight="bold", pad=1.4)
        r = thresholds.loc[city]
        if r.threshold_status == "no_positive_gain":
            ax.set_title("No positive gain for women", loc="right", color="#702d8a", fontsize=6.4, pad=2.5)
            continue
        onset, l90 = 1000 * r.female_onset, 1000 * r.female_L90
        for j, value in enumerate([onset, l90]):
            coincident = onset == l90
            if value <= 12:
                start, end = ((0, .5) if j == 0 else (.5, 1)) if coincident else (0, 1)
                ax.axvline(value, start, end, color=SEX["female"], lw=.6, ls=":" if j == 0 else "-")
            if j == 1 and coincident:
                continue
            at = min(value, 12)
            color = SEX["female"]
            marker = ">" if value > 12 else ("o" if j == 0 else "|")
            ax.plot(at, .025, marker, color=color, ms=2.5, transform=ax.get_xaxis_transform(), clip_on=False)
            align = "right" if at > 2.95 else "left"
            if city in ["hong_kong", "sydney"] and j == 1:
                align = "left"
            if city == "auckland" and j == 1:
                align = "right"
            label = f"{value:.2f}" + (" →" if value > 12 else "")
            ax.annotate(label, (at, label_heights[index][j]), xycoords=ax.get_xaxis_transform(),
                        xytext=(2 if align == "left" else -2, 0), textcoords="offset points",
                        ha=align, va="top", color="#762f92", fontsize=7.8)
    text_at(fig, 30, 1141, "Nighttime PA relative to the local\nmedian (×, log scale)",
            ha="center", rotation=90, fontsize=8.5)
    handles = [Line2D([], [], color=c, lw=1.2) for c in SEX.values()]
    handles += [Patch(color=SEX["female"], alpha=.17, lw=0)]
    handles += [Line2D([], [], color=SEX["female"], ls=ls, lw=.7) for ls in [":", "-"]]
    fig.legend(handles, ["Women", "Men", "95% CI", "Women’s onset", "Women’s L90"],
               loc="center right", bbox_to_anchor=(.971, .029), ncol=5, fontsize=6.7,
               handlelength=1.7, columnspacing=1.0)
    text_at(fig, 563, 1480, RAD, ha="center", fontsize=9)


def main():
    fig = plt.figure(figsize=(8.5, 8.5 * 1500 / 1133))
    top_panels(fig)
    context_panels(fig)
    city_panels(fig)
    for label, x, y in [("a", 36, 30), ("b", 516, 30), ("c", 58, 515),
                        ("d", 624, 515), ("e", 58, 820)]:
        text_at(fig, x, y, label, fontsize=14.5, weight="bold")
    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    for suffix in ["png", "pdf", "svg"]:
        metadata = {"Date": None} if suffix == "svg" else None
        fig.savefig(output / f"fig3.{suffix}", dpi=450, metadata=metadata)
    plt.close(fig)


if __name__ == "__main__":
    main()
