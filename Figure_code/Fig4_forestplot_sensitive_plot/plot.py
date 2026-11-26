"""Assemble the published Fig. 4 values into the a, d, e, f, g layout."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
import pandas as pd

ROOT = Path(__file__).resolve().parent
COLOURS = {'female': '#8E3B9A', 'male': '#E17E3C'}
OFFSETS = {'female': -.17, 'male': .17}
SCOPES = {'all25': 'All 25 cities', 'recreational': 'Recreational',
          'residential': 'Residential', 'commercial': 'Commercial'}
GROUPS = {
    'Population & socioeconomic conditions': 'Population &\nsocioeconomic\nconditions',
    'City lighting structure': 'City lighting\nstructure',
    'City social context': 'City social\ncontext',
    'Climate & temporal effects': 'Climate &\ntemporal effects',
}
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 13, 'pdf.fonttype': 42,
                    'svg.fonttype': 'none', 'axes.spines.top': False,
                    'axes.spines.left': False, 'axes.spines.right': False,
                    'axes.edgecolor': '#5b5b5b', 'axes.linewidth': .65,
                    'mathtext.default': 'regular', 'savefig.facecolor': 'white'})


def load(name):
    return pd.read_csv(ROOT / 'data' / name, encoding='utf-8-sig')


def letter(fig, mark, x, y):
    fig.text(x, y, mark, size=29, weight='bold', va='top')


def errorbar(ax, value, lo, hi, y, colour, ms=5):
    assert lo <= value <= hi
    ax.errorbar(value, y, xerr=[[value-lo], [hi-value]], capsize=2.2,
                elinewidth=1.15, marker='o', ms=ms, color=colour)


def coefficient_panel(fig, coefficients):
    df = pd.read_csv(coefficients)
    female = df[df.sex.eq('female')].sort_values('plot_order')
    ordered = pd.concat([female[female.term.eq('z_road_mean')], female[~female.term.eq('z_road_mean')]])
    terms = list(ordered.term)
    assert len(df) == 62 and len(set(terms)) == 31 and terms[-1] == 'is_weekend'
    start, bottom, right, top = .018, .158, .743, .935
    table = fig.add_axes([start, bottom, right-start, top-bottom], frameon=False)
    forest = fig.add_axes([.300, bottom, .182, top-bottom])
    forest.patch.set_visible(False)
    for ax in (table, forest):
        ax.set_ylim(30.5, -.5)
        ax.set_yticks([])
    table.set_xticks([])
    table.set_xlim(start, right)
    forest.set_xlim(-1.95, 1.78)
    forest.set_xticks([-1, 0, 1])
    forest.axvline(0, color='#8c8c8c', lw=.8, ls='--', zorder=0)
    forest.tick_params(axis='x', labelsize=13, pad=3)
    forest.set_xlabel('Log coefficient (95% CI)', size=14, labelpad=4)
    letter(fig, 'a', start, .983)
    for i, row in enumerate(ordered.itertuples(index=False)):
        table.text(.111, i, row.label, va='center', size=13.8)
    blocks = ordered.block.tolist()
    for block in dict.fromkeys(blocks):
        rows = [i for i, value in enumerate(blocks) if value == block]
        group = GROUPS.get(block, block.replace('CPTED ', 'CPTED\n') if block != 'CPTED concealment' else block)
        table.text(.021, (rows[0]+rows[-1])/2, group, weight='bold', va='center', size=12.8, linespacing=1.0)
        if rows[-1] != 30:
            table.axhline(rows[-1]+.5, color='#c7c7c7', lw=.6, zorder=0)
    for sex, anchor in [('female', .504), ('male', .634)]:
        fig.text(anchor, .949, f'{sex.title()} beta [95% CI]', color=COLOURS[sex], weight='bold', size=14)
        records = df[df.sex.eq(sex)].set_index('term').loc[terms]
        for i, row in enumerate(records.itertuples()):
            errorbar(forest, row.estimate, row.lo_95_normal, row.hi_95_normal, i+OFFSETS[sex], COLOURS[sex], 5.1)
            stars = '***' if row.p < .001 else '**' if row.p < .01 else '*' if row.p < .05 else ''
            table.text(anchor-.0035, i, stars, va='center', ha='right', size=12.7)
            table.text(anchor, i, f'{row.estimate:+.3f} [{row.lo_95_normal:+.3f}, {row.hi_95_normal:+.3f}]',
                       va='center', size=13)
    fig.text(.021, .097, 'Model-based 95% CIs and P values.  * P < 0.05   ** P < 0.01   *** P < 0.001',
             size=12.1, va='center')
    legend = [Line2D([], [], color=colour, marker='o', ls='none', label=sex.title())
              for sex, colour in COLOURS.items()]
    fig.legend(handles=legend, loc='center', bbox_to_anchor=(.63, .097), ncol=2,
               handletextpad=.2, frameon=False, fontsize=12.5)
    return table, forest


def side_axis(fig, mark, title, rect, title_y, labels, reference):
    letter(fig, mark, .764, title_y+.007)
    fig.text(.795, title_y, title, weight='bold', va='top', size=14.5)
    ax = fig.add_axes(rect)
    ax.set_ylim(len(labels)-.55, -.45)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, size=11.7)
    ax.tick_params(axis='y', length=0, pad=5)
    ax.tick_params(axis='x', length=3, labelsize=11.5, pad=3)
    ax.axvline(reference, color='#9c9c9c', lw=.7, ls='--', zorder=0)
    return ax


def ratio_scale(ax, limits, ticks):
    ax.set_xscale('log')
    ax.set_xlim(*limits)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:g}' for t in ticks])
    ax.xaxis.set_minor_locator(NullLocator())


def side_panels(fig):
    labels = list(SCOPES.values())
    groups = [(scope, period) for scope in SCOPES for period in ('night', 'day')]
    dist_labels = [f'{SCOPES[scope]} · {period}' for scope, period in groups]
    ax_d = side_axis(fig, 'd', 'Night / day lighting gain', [.838, .798, .145, .142], .974, labels, 1)
    ax_e = side_axis(fig, 'e', 'Clock versus solar timing', [.838, .583, .145, .142], .758, labels, 1)
    ax_f = side_axis(fig, 'f', 'Change in fitted light-response shape', [.850, .282, .133, .211], .542, dist_labels, 0)
    ax_g = side_axis(fig, 'g', 'Low / high safety-tercile lighting gain', [.838, .056, .129, .135], .235, labels, 1)
    ratio_scale(ax_d, (.60, 7), [.75, 1, 2, 4, 6])
    ax_d.set_xlabel('Gain ratio (95% paired FRW)', size=11.8)
    ax_e.set_xlim(.65, 6)
    ax_e.set_xticks(range(1,7))
    ax_e.set_xlabel('Gain ratio (hollow: clock, filled: solar)', size=11.8)
    ax_f.set_xlim(-.05, 5.15)
    ax_f.set_xticks(range(6))
    ax_f.set_xlabel('Curve distance from covered baseline (%)', size=11.8)
    ratio_scale(ax_g, (.1, 3.4), [.1, .25, .5, 1, 2])
    ax_g.set_xlabel('Night gain ratio (95% paired FRW)', size=11.8)
    handles = [Line2D([], [], color='#45434a', marker=marker, mfc=face, ls='none', label=label)
               for marker, face, label in [('o', '#45434a', '+ Safer'), ('D', 'white', '+ Safer and livelier')]]
    ax_f.legend(handles=handles, loc='lower left', bbox_to_anchor=(0,1.005), frameon=False,
                ncol=2, fontsize=10.8, handletextpad=.3, columnspacing=1, borderaxespad=0)
    gain = load('night_day_gain.csv').set_index(['scope', 'sex'])
    solar = load('clock_solar.csv').set_index(['scope', 'sex'])
    perception = load('perception_adjustment.csv').set_index(['scope', 'period', 'sex'])
    safe = load('safety_terciles.csv').set_index(['scope', 'contrast'])
    assert len(perception) == 16 and perception.index.is_unique
    for y, scope in enumerate(SCOPES):
        for sex, offset in OFFSETS.items():
            c = COLOURS[sex]
            row = gain.loc[(scope, sex)]
            errorbar(ax_d, row.night_over_day_gain_ratio, row.lo, row.hi, y+offset, c)
            row = solar.loc[(scope, sex)]
            ax_e.plot([row.clock, row.solar], [y+offset]*2, lw=1.2, color=c, alpha=.57)
            ax_e.plot(row.clock, y+offset, 'o', color=c, mfc='white', ms=5)
            ax_e.plot(row.solar, y+offset, 'o', color=c, ms=5)
            row = safe.loc[(scope, sex+' low/high')]
            errorbar(ax_g, row.ratio, row.lo, row.hi, y+offset, c)
    for y, (scope, period) in enumerate(groups):
        for sex, offset in OFFSETS.items():
            c = COLOURS[sex]
            row = perception.loc[(scope, period, sex)]
            safer = 100*row.reldiff_pp_vs_sub
            joint = 100*row.reldiff_pp2_vs_sub
            ax_f.plot([safer, joint], [y+offset]*2, color=c, lw=1.0, alpha=.60)
            ax_f.plot(safer, y+offset, 'o', color=c, ms=4.8)
            ax_f.plot(joint, y+offset, 'D', color=c, mfc='white', ms=4.8)
    return ax_d, ax_e, ax_f, ax_g


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refit', action='store_true', help='Use panel-a estimates from model/fit.R.')
    args = parser.parse_args()
    coefficients = ROOT / ('output/refit_coefficients.csv' if args.refit else 'data/coefficients.csv')
    fig = plt.figure(figsize=(21, 12.35))
    coefficient_panel(fig, coefficients)
    side_panels(fig)
    output = ROOT / 'output'
    output.mkdir(exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(output / f'fig4.{ext}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
