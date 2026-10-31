#!/usr/bin/env python3
"""
12_main_panels.py  (v19)

Changes (relative to v18):
  - fmt_rho now uses dynamic decimal places:
      |rho| < 0.005   → display with 4 decimal places (e.g. "-0.0003"), preserving direction
      |rho| exactly 0 → "0.00"
      otherwise       → 2 decimal places (e.g. "+0.11")
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import MaxNLocator


# Global constants
MAX_CTRAC_DISPLAY = 6
JITTER_A = 0.30


# Nature style
def set_nature_style():
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 11.5,
        'axes.linewidth': 0.9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.labelsize': 9.5,
        'ytick.labelsize': 10,
        'xtick.major.width': 0.9,
        'ytick.major.width': 0.9,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'legend.fontsize': 9,
        'legend.frameon': False,
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.08,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'pdf.compression': 6,
    })


COLORS_CTRAC = {'C3': '#A6CEE3', 'C4': '#4A90C2', 'C5+': '#1F4E79'}
COLOR_LOOP = '#7A7A7A'
COLOR_CDEN = '#2E7D32'
COLOR_FIT  = '#E53935'


# Utilities
def p_to_star(p):
    if np.isnan(p): return 'ns'
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'


def fmt_p(p):
    return f'{p:.1e}' if p < 1e-3 else f'{p:.3f}'


def fmt_rho(rho, ndigits=2, eps_small=0.005, ndigits_small=4):
    """
    Format rho with dynamic decimal places:
      - |rho| < 10^-4          → "0.00" (true zero)
      - |rho| < eps_small      → ndigits_small decimal places, preserving sign
                                 e.g. -0.000288 → "-0.0003"
      - otherwise              → ndigits decimal places, preserving sign
                                 e.g. +0.109826 → "+0.11"
    """
    if rho is None or (isinstance(rho, float) and np.isnan(rho)):
        return 'nan'
    if abs(rho) < 10 ** (-ndigits_small):
        return f'{0.0:.{ndigits}f}'
    if abs(rho) < eps_small:
        return f'{rho:+.{ndigits_small}f}'
    return f'{rho:+.{ndigits}f}'


def filter_iqr(x, y, iqr_k=2.0):
    if iqr_k is None or iqr_k <= 0:
        return x, y
    mask = np.ones(len(x), dtype=bool)
    for arr in (x, y):
        q1, q3 = np.percentile(arr, [25, 75])
        iqr = q3 - q1
        if iqr <= 0:
            continue
        mask &= (arr >= q1 - iqr_k * iqr) & (arr <= q3 + iqr_k * iqr)
    return x[mask], y[mask]


def spearman_with_ci(x, y, alpha=0.05):
    n = len(x)
    if n < 5:
        return np.nan, np.nan, np.nan, np.nan
    rho, p = stats.spearmanr(x, y)
    if abs(rho) >= 1.0:
        return rho, p, rho, rho
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3) if n > 3 else np.nan
    z_crit = stats.norm.ppf(1 - alpha / 2)
    ci_low = np.tanh(z - z_crit * se)
    ci_high = np.tanh(z + z_crit * se)
    return rho, p, ci_low, ci_high


# Common sample mask
def get_common_mask(df, panels, iqr_k=2.0):
    mask = df['pHT'].notna().values.copy()
    for _, col, _ in panels:
        mask &= df[col].notna().values
    n_base = int(mask.sum())

    if iqr_k <= 0 or n_base == 0:
        return mask, n_base, n_base

    for _, col, _ in panels:
        idx = np.where(mask)[0]
        if len(idx) == 0:
            break
        x = df.iloc[idx][col].values.astype(float)
        y = df.iloc[idx]['pHT'].values.astype(float)
        keep = np.ones(len(idx), dtype=bool)
        for arr in (x, y):
            q1, q3 = np.percentile(arr, [25, 75])
            iqr = q3 - q1
            if iqr <= 0:
                continue
            keep &= (arr >= q1 - iqr_k * iqr) & (arr <= q3 + iqr_k * iqr)
        mask[idx[~keep]] = False

    return mask, n_base, int(mask.sum())


# Combined box in the upper right
def add_combined_box_side_by_side(ax,
                                   colors, labels,
                                   rho, p, n,
                                   anchor=(0.97, 0.97),
                                   fontsize=8.5,
                                   box_w=0.50,
                                   pad_x=0.028,
                                   pad_y=0.024,
                                   row_h=0.052,
                                   dot_offset_x=0.018,
                                   text_offset_x=0.048,
                                   stat_offset_ratio=0.52):
    n_rows = max(len(labels), 3)
    box_h = pad_y * 2 + row_h * n_rows
    x1 = anchor[0] - box_w
    y1 = anchor[1] - box_h

    box = FancyBboxPatch(
        (x1, y1), box_w, box_h,
        boxstyle='round,pad=0.005,rounding_size=0.012',
        transform=ax.transAxes,
        facecolor='white', edgecolor='#CCCCCC',
        linewidth=0.8, zorder=9)
    ax.add_patch(box)

    for i, (color, label) in enumerate(zip(colors, labels)):
        y = anchor[1] - pad_y - row_h * (i + 0.5)
        x_dot = x1 + pad_x + dot_offset_x
        x_text = x1 + pad_x + text_offset_x
        ax.plot(x_dot, y, marker='o',
                markersize=6.5, color=color,
                markeredgecolor='white', markeredgewidth=0.7,
                transform=ax.transAxes, zorder=10, clip_on=False)
        ax.text(x_text, y, label,
                transform=ax.transAxes,
                ha='left', va='center',
                fontsize=fontsize, zorder=10)

    stat_x = x1 + box_w * stat_offset_ratio
    stats_text = [
        f'n = {n}',
        f'\u03c1 = {fmt_rho(rho)}',
        f'P = {fmt_p(p)}',
    ]
    for i, text in enumerate(stats_text):
        y = anchor[1] - pad_y - row_h * (i + 0.5)
        ax.text(stat_x, y, text,
                transform=ax.transAxes,
                ha='left', va='center',
                fontsize=fontsize, zorder=10)


# Panel a
def plot_panel_a(df, ax):
    sub = df.dropna(subset=['ctracNum', 'pHT']).copy()

    def classify(x):
        if x <= 3: return 'C3'
        elif x == 4: return 'C4'
        else: return 'C5+'

    sub['ctrac_class'] = sub['ctracNum'].apply(classify)
    sub['ctrac_disp'] = sub['ctracNum'].clip(upper=MAX_CTRAC_DISPLAY)
    n_collapsed = int((sub['ctracNum'] > MAX_CTRAC_DISPLAY).sum())

    rng = np.random.default_rng(42)
    for cls in ['C3', 'C4', 'C5+']:
        s = sub[sub['ctrac_class'] == cls]
        if len(s) == 0:
            continue
        x = s['ctrac_disp'].values.astype(float)
        xj = x + rng.uniform(-JITTER_A, JITTER_A, len(x))
        ax.scatter(xj, s['pHT'],
                   s=24, alpha=0.72, color=COLORS_CTRAC[cls],
                   edgecolor='white', linewidth=0.45,
                   zorder=3)

    x_all = sub['ctrac_disp'].values.astype(float)
    y_all = sub['pHT'].values
    z = np.polyfit(x_all, y_all, 1)
    xs = np.linspace(x_all.min(), x_all.max(), 100)
    ax.plot(xs, np.polyval(z, xs),
            color=COLOR_FIT, lw=1.6, ls='--',
            zorder=4, alpha=0.9)

    rho, p = stats.spearmanr(sub['ctracNum'].values, sub['pHT'].values)

    add_combined_box_side_by_side(
        ax,
        colors=[COLORS_CTRAC['C3'],
                COLORS_CTRAC['C4'],
                COLORS_CTRAC['C5+']],
        labels=['C3', 'C4', 'C5+'],
        rho=rho, p=p, n=len(sub),
        anchor=(0.97, 0.97),
        fontsize=8.5, box_w=0.50)

    x_data_max = int(np.ceil(sub['ctrac_disp'].max()))
    ax.set_xlim(2.5 - JITTER_A, x_data_max + 0.5 + JITTER_A)
    ax.set_xticks(np.arange(3, x_data_max + 1))

    if n_collapsed > 0:
        ax.text(0.02, 0.02,
                f'C-tract >{MAX_CTRAC_DISPLAY} collapsed (n={n_collapsed})',
                transform=ax.transAxes,
                fontsize=7, color='#888888',
                ha='left', va='bottom')

    y_min, y_max = y_all.min(), y_all.max()
    y_range = y_max - y_min
    ax.set_ylim(y_min - y_range * 0.03, y_max + y_range * 0.22)

    ax.set_xlabel('C-tract length (stem)', fontsize=11)
    ax.set_ylabel('pHT', fontsize=11)
    ax.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title('C-tract length vs pHT',
                 fontsize=11.5, fontweight='bold', loc='left', pad=8)

    return rho, p, len(sub)


# Forest plot
def draw_forest(df, panels, ax, color,
                title='', xlim=None,
                null_band=0.10,
                show_ci=True,
                iqr_k=2.0,
                common_mask=None,
                common_n_base=None):
    rows = []
    for letter, col, label in panels:
        if common_mask is not None:
            sub = df[common_mask]
            sub = sub.dropna(subset=[col, 'pHT'])
            x_f = sub[col].values.astype(float)
            y_f = sub['pHT'].values.astype(float)
            n_raw = common_n_base if common_n_base is not None else len(sub)
        else:
            sub = df.dropna(subset=[col, 'pHT'])
            if len(sub) < 5:
                continue
            x_raw = sub[col].values.astype(float)
            y_raw = sub['pHT'].values.astype(float)
            n_raw = len(x_raw)
            x_f, y_f = filter_iqr(x_raw, y_raw, iqr_k=iqr_k)

        if len(x_f) < 5:
            continue
        rho, p, ci_low, ci_high = spearman_with_ci(x_f, y_f)
        rows.append({
            'letter': letter, 'feature': col, 'label': label,
            'n': len(x_f), 'n_raw': n_raw,
            'rho': rho, 'p': p,
            'ci_low': ci_low, 'ci_high': ci_high,
            'star': p_to_star(p),
        })

    stats_df = pd.DataFrame(rows)
    if stats_df.empty:
        return stats_df

    stats_df = stats_df.iloc[::-1].reset_index(drop=True)
    y_pos = np.arange(len(stats_df))

    if xlim is None:
        all_vals = np.concatenate([
            stats_df['ci_low'].values,
            stats_df['ci_high'].values,
            [-0.05, 0.05]])
        xmin = min(all_vals) - 0.05
        xmax = max(all_vals) + 0.05
    else:
        xmin, xmax = xlim
    x_span = xmax - xmin

    for i in y_pos:
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color='#F7F7F7', zorder=0)

    ax.axvspan(-null_band, null_band,
               color='#EFEFEF', alpha=0.7, zorder=0)
    ax.axvline(0, color='#333333', lw=1.2, ls='-', alpha=0.7, zorder=2)

    if show_ci:
        for i, row in stats_df.iterrows():
            ax.plot([row['ci_low'], row['ci_high']], [i, i],
                    color=color, lw=2.0, alpha=0.45,
                    solid_capstyle='round', zorder=3)

    ax.scatter(stats_df['rho'], y_pos,
               s=75, color=color,
               edgecolor='white', linewidth=1.2, zorder=4)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(stats_df['label'], fontsize=9.5)
    ax.set_ylim(-0.6, len(stats_df) - 0.4)

    x_text_left = xmax + x_span * 0.08
    x_text_right = xmax + x_span * 0.42
    ax.set_xlim(xmin, xmax + x_span * 0.55)

    ax.set_xlabel('Spearman \u03c1', fontsize=11, labelpad=6)
    ax.grid(alpha=0.18, ls='--', axis='x', lw=0.6, zorder=1)
    ax.set_axisbelow(False)

    for i, row in stats_df.iterrows():
        # Use fmt_rho (dynamic decimal places)
        ax.text(x_text_left, i, fmt_rho(row['rho']),
                ha='left', va='center',
                fontsize=9.5, color='#333333')
        star_color = '#333333' if row['star'] != 'ns' else '#999999'
        ax.text(x_text_right, i, row['star'],
                ha='left', va='center',
                fontsize=10.5, fontweight='bold', color=star_color)

    if title:
        ax.set_title(title, fontsize=11.5,
                     fontweight='bold', loc='left', pad=10)

    return stats_df


LOOP_PANELS = [
    ('a', 'maxLoop',           'Max loop length'),
    ('b', 'minLoop',           'Min loop length'),
    ('c', 'allLoop_len',       'All loops sum'),
    ('d', 'loop2_len',         'Middle loop length'),
    ('e', 'side_loop_sum',     'Side loops sum'),
    ('f', 'side_loop_max_len', 'Longest side loop'),
    ('g', 'side_loop_min_len', 'Shortest side loop'),
]

CDEN_PANELS = [
    ('a', 'all_C_density',           'C density (iM)'),
    ('b', 'allLoop_C_density',       'C density (all loops)'),
    ('c', 'loop2_C_density',         'C density (middle loop)'),
    ('d', 'side_loop_C_density',     'C density (side loops)'),
    ('e', 'max_side_loop_C_density', 'C density (longest side)'),
    ('f', 'min_side_loop_C_density', 'C density (shortest side)'),
]


def plot_panel_b_forest(df, ax, iqr_k=2.0,
                        common_mask=None, common_n_base=None):
    return draw_forest(
        df, LOOP_PANELS, ax, COLOR_LOOP,
        title='Loop length metrics vs pHT',
        xlim=(-0.25, 0.35),
        iqr_k=iqr_k,
        common_mask=common_mask,
        common_n_base=common_n_base)


def plot_panel_c_forest(df, ax, iqr_k=0.0):
    return draw_forest(
        df, CDEN_PANELS, ax, COLOR_CDEN,
        title='C density metrics vs pHT',
        xlim=(-0.05, 0.60),
        iqr_k=iqr_k)


# Single-feature scatter plots
def draw_single_scatter_panel(df, col, label, ax, color,
                               letter=None,
                               integer_x=False,
                               iqr_k=2.0,
                               common_mask=None):
    if common_mask is not None:
        sub = df[common_mask].dropna(subset=[col, 'pHT']).copy()
        x = sub[col].values.astype(float)
        y = sub['pHT'].values.astype(float)
    else:
        sub = df.dropna(subset=[col, 'pHT']).copy()
        if len(sub) < 5:
            ax.text(0.5, 0.5, 'insufficient data',
                    transform=ax.transAxes, ha='center', va='center',
                    color='#999999')
            return None
        x_raw = sub[col].values.astype(float)
        y_raw = sub['pHT'].values.astype(float)
        x, y = filter_iqr(x_raw, y_raw, iqr_k=iqr_k)

    if len(x) < 5:
        ax.text(0.5, 0.5, 'insufficient after filter',
                transform=ax.transAxes, ha='center', va='center',
                color='#999999')
        return None

    rho, p = stats.spearmanr(x, y)

    rng = np.random.default_rng(42)
    span = x.max() - x.min() if x.max() > x.min() else 1.0
    x_plot = x + rng.uniform(-span * 0.012, span * 0.012, len(x))

    ax.scatter(x_plot, y,
               s=20, alpha=0.65, color=color,
               edgecolor='white', linewidth=0.4, zorder=3)

    if x.max() > x.min():
        z = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, np.polyval(z, xs),
                color=COLOR_FIT, lw=1.3, ls='--',
                alpha=0.85, zorder=4)

    box_text = f'n = {len(x)}\n\u03c1 = {fmt_rho(rho)}\nP = {fmt_p(p)}'
    ax.text(0.97, 0.97, box_text,
            transform=ax.transAxes, ha='right', va='top',
            fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.35',
                      facecolor='white', edgecolor='#CCCCCC',
                      linewidth=0.8),
            zorder=10)

    ax.set_xlabel(label, fontsize=10)
    ax.set_ylabel('pHT', fontsize=10)
    ax.grid(alpha=0.18, ls='--', lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    if integer_x:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

    if letter:
        ax.text(-0.22, 1.05, letter,
                transform=ax.transAxes,
                fontsize=14, fontweight='bold',
                va='bottom', ha='left')

    return rho, p, len(x)


def plot_scatter_grid(df, panels, outdir, prefix, color,
                      ncols=4, integer_x=False, dpi=300,
                      iqr_k=2.0, common_mask=None):
    n_panels = len(panels)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.8 * ncols, 3.6 * nrows),
        squeeze=False)

    results = []
    for i, (letter, col, label) in enumerate(panels):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        res = draw_single_scatter_panel(
            df, col, label, ax, color,
            letter=letter, integer_x=integer_x,
            iqr_k=iqr_k, common_mask=common_mask)
        if res is not None:
            results.append({
                'letter': letter, 'feature': col, 'label': label,
                'n': res[2], 'rho': res[0], 'p': res[1],
                'star': p_to_star(res[1]),
            })

    for j in range(n_panels, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis('off')

    fig.tight_layout()
    fig.savefig(outdir / f'{prefix}.png', dpi=dpi, bbox_inches='tight')
    fig.savefig(outdir / f'{prefix}.pdf', dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    return pd.DataFrame(results)


# main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--iqr-k', type=float, default=2.0,
                    help='IQR filter coefficient for Panel b (default 2.0).')
    ap.add_argument('--common-n', dest='common_n',
                    action='store_true', default=True,
                    help='Use the same set of samples for all features in Panel b (default on).')
    ap.add_argument('--no-common-n', dest='common_n',
                    action='store_false',
                    help='Disable common samples and revert to independent filtering per feature.')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep='\t')
    df = df[df['ctracNum'].notna()].reset_index(drop=True)
    print(f'Loaded {len(df)} valid iM sequences')
    print(f'pHT samples: {df["pHT"].notna().sum()}')
    print(f'Panel b IQR k = {args.iqr_k} '
          f'({"on" if args.iqr_k > 0 else "off"})')
    print(f'Panel b common samples: {"on" if args.common_n else "off"}')
    print(f'Panel c filtering: off')

    set_nature_style()

    # Panel a
    fig_a, ax_a = plt.subplots(figsize=(5.2, 5))
    rho_a, p_a, n_a = plot_panel_a(df, ax_a)
    fig_a.savefig(outdir / 'fig2_main_a.png', dpi=args.dpi, bbox_inches='tight')
    fig_a.savefig(outdir / 'fig2_main_a.pdf', dpi=args.dpi, bbox_inches='tight')
    plt.close(fig_a)
    print(f'\nsaved: fig2_main_a.png')
    print(f'  a: n={n_a}, \u03c1={rho_a:+.3f}, P={p_a:.2e}')

    # Panel b common sample mask
    if args.common_n:
        common_mask_b, n_base_b, n_final_b = get_common_mask(
            df, LOOP_PANELS, iqr_k=args.iqr_k)
        print(f'\nPanel b common samples: '
              f'common non-null {n_base_b} → IQR intersection {n_final_b}')
    else:
        common_mask_b, n_base_b = None, None

    # Panel b: forest plot
    fig_b, ax_b = plt.subplots(figsize=(7.5, 5))
    stats_b = plot_panel_b_forest(
        df, ax_b, iqr_k=args.iqr_k,
        common_mask=common_mask_b, common_n_base=n_base_b)
    fig_b.savefig(outdir / 'fig2_main_b.png', dpi=args.dpi, bbox_inches='tight')
    fig_b.savefig(outdir / 'fig2_main_b.pdf', dpi=args.dpi, bbox_inches='tight')
    plt.close(fig_b)
    print(f'\nsaved: fig2_main_b.png')
    print(stats_b[['label', 'n', 'n_raw', 'rho', 'p', 'star']]
          .to_string(index=False))
    stats_b.to_csv(outdir / 'fig2_main_b_stats.csv', index=False)

    stats_b_scatter = plot_scatter_grid(
        df, LOOP_PANELS, outdir,
        prefix='fig2_main_b_scatter',
        color=COLOR_LOOP,
        ncols=4, integer_x=True,
        dpi=args.dpi, iqr_k=args.iqr_k,
        common_mask=common_mask_b)
    print(f'\nsaved: fig2_main_b_scatter.png (7 panels)')
    print(stats_b_scatter[['label', 'n', 'rho', 'p', 'star']]
          .to_string(index=False))
    stats_b_scatter.to_csv(
        outdir / 'fig2_main_b_scatter_stats.csv', index=False)

    # Panel c: forest plot (no filtering)
    fig_c, ax_c = plt.subplots(figsize=(7.5, 4.5))
    stats_c = plot_panel_c_forest(df, ax_c, iqr_k=0.0)
    fig_c.savefig(outdir / 'fig2_main_c.png', dpi=args.dpi, bbox_inches='tight')
    fig_c.savefig(outdir / 'fig2_main_c.pdf', dpi=args.dpi, bbox_inches='tight')
    plt.close(fig_c)
    print(f'\nsaved: fig2_main_c.png')
    print(stats_c[['label', 'n', 'n_raw', 'rho', 'p', 'star']]
          .to_string(index=False))
    stats_c.to_csv(outdir / 'fig2_main_c_stats.csv', index=False)

    stats_c_scatter = plot_scatter_grid(
        df, CDEN_PANELS, outdir,
        prefix='fig2_main_c_scatter',
        color=COLOR_CDEN,
        ncols=4, integer_x=False,
        dpi=args.dpi, iqr_k=0.0)
    print(f'\nsaved: fig2_main_c_scatter.png (6 panels)')
    print(stats_c_scatter[['label', 'n', 'rho', 'p', 'star']]
          .to_string(index=False))
    stats_c_scatter.to_csv(
        outdir / 'fig2_main_c_scatter_stats.csv', index=False)

    # Three panels in one row
    fig_all = plt.figure(figsize=(20, 5.5))
    gs = fig_all.add_gridspec(
        1, 3,
        width_ratios=[1.0, 1.35, 1.35],
        wspace=0.55,
        left=0.04, right=0.98,
        top=0.88, bottom=0.13)

    ax_a = fig_all.add_subplot(gs[0, 0])
    ax_b = fig_all.add_subplot(gs[0, 1])
    ax_c = fig_all.add_subplot(gs[0, 2])

    plot_panel_a(df, ax_a)
    plot_panel_b_forest(
        df, ax_b, iqr_k=args.iqr_k,
        common_mask=common_mask_b, common_n_base=n_base_b)
    plot_panel_c_forest(df, ax_c, iqr_k=0.0)

    for ax, letter in zip([ax_a, ax_b, ax_c], ['E', 'F', 'G']):
        ax.text(-0.16, 1.06, letter,
                transform=ax.transAxes,
                fontsize=20, fontweight='bold',
                va='bottom', ha='left')

    fig_all.savefig(outdir / 'fig2_main_combined.png',
                    dpi=args.dpi, bbox_inches='tight')
    fig_all.savefig(outdir / 'fig2_main_combined.pdf',
                    dpi=args.dpi, bbox_inches='tight')
    plt.close(fig_all)
    print(f'\nsaved: fig2_main_combined.png')
    print(f'\nDone, results in {outdir}')


if __name__ == '__main__':
    main()
