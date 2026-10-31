#!/usr/bin/env python3
"""
0333_tss_enrichment.py  (v11 - 主图去掉 CI，只画曲线 + 填充)

- combined_TSS_facet（主图）：单样本曲线，无 CI
- combined_TSS_species_mean_facet（附图）：均值 ± SEM + 统计框
"""

import os
import gzip
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter1d
from scipy.stats import ks_2samp

# Nature style
def set_nature_style():
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 11,
        'axes.labelsize': 11.5,
        'axes.titlesize': 12,
        'axes.linewidth': 0.9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'xtick.major.width': 0.9,
        'ytick.major.width': 0.9,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'legend.fontsize': 9.5,
        'legend.frameon': False,
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'pdf.compression': 6,
    })

# Colors
COLORS_HUMAN = {
    'HEK293T_iM_highConf_intersect': '#08519C',
    'WDLPS_iM_highConf_intersect':   '#9ECAE1',
}
COLORS_RICE = {
    'IP_vs_IgG_intersect':   '#A1D99B',
    'IP_vs_Input_intersect': '#238B45',
}
COLOR_HUMAN_MEAN = '#08519C'
COLOR_RICE_MEAN  = '#238B45'

HUMAN_DISPLAY = {
    'HEK293T_iM_highConf_intersect': 'Human CUT&Tag HEK293T',
    'WDLPS_iM_highConf_intersect':   'Human CUT&Tag WDLPS',
}
RICE_DISPLAY = {
    'IP_vs_IgG_intersect':   'Rice iM-IP (CK pH7.0)',
    'IP_vs_Input_intersect': 'Rice iM-IP (seedling)',
}
HUMAN_GROUPS = [
    ('Human CUT&Tag', [
        'HEK293T_iM_highConf_intersect',
        'WDLPS_iM_highConf_intersect',
     ], COLOR_HUMAN_MEAN),
]
RICE_GROUPS = [
    ('Rice iM-IP', [
        'IP_vs_IgG_intersect',
        'IP_vs_Input_intersect',
     ], COLOR_RICE_MEAN),
]

# IO
def open_maybe_gz(path):
    return gzip.open(str(path), 'rt') if str(path).endswith('.gz') else open(path)

def read_peak_bed(path, min_qvalue=None):
    rows = []
    with open_maybe_gz(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            try:
                chrom = parts[0]; start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            if min_qvalue is not None and len(parts) >= 9:
                try:
                    qval = float(parts[8])
                except ValueError:
                    continue
                if qval < min_qvalue:
                    continue
            rows.append((chrom, start, end))
    if not rows:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])

def read_genes_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            try:
                rows.append((parts[0], int(parts[1]), int(parts[2]), parts[5]))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'strand'])

def read_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                sizes[p[0]] = int(p[1])
    return sizes

def merged_length(df):
    if df.empty:
        return 0
    total = 0
    for _, grp in df.groupby('chrom', sort=False):
        ivs = sorted(zip(grp['start'].tolist(), grp['end'].tolist()))
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                total += ce - cs
                cs, ce = s, e
        total += ce - cs
    return total

def strip_suffix(name):
    for suf in ('.narrowPeak.gz', '.narrowPeak', '.bed.gz', '.bed'):
        if name.endswith(suf):
            return name[:-len(suf)]
    return name

def build_peak_lookup(peak_df):
    out = {}
    for chrom, grp in peak_df.groupby('chrom', sort=False):
        s = grp['start'].to_numpy(dtype=np.int64)
        e = grp['end'].to_numpy(dtype=np.int64)
        order = np.argsort(s)
        out[chrom] = (s[order], e[order])
    return out

def accumulate_tss(genes_df, lookup, window):
    length = 2 * window + 1
    diff = np.zeros(length + 1, dtype=np.int64)
    n_genes = 0
    for chrom, gs, ge, strand in zip(genes_df['chrom'],
                                     genes_df['start'],
                                     genes_df['end'],
                                     genes_df['strand']):
        if strand == '+':
            tss = int(gs)
        elif strand == '-':
            tss = int(ge)
        else:
            continue
        n_genes += 1
        if chrom not in lookup:
            continue
        s_arr, e_arr = lookup[chrom]
        qs = tss - window
        qe = tss + window + 1
        hi = np.searchsorted(s_arr, qe, side='left')
        if hi == 0:
            continue
        cand_s = s_arr[:hi]
        cand_e = e_arr[:hi]
        mask = cand_e > qs
        if not mask.any():
            continue
        cand_s = cand_s[mask]
        cand_e = cand_e[mask]
        if strand == '+':
            rel_s = cand_s - tss
            rel_e = cand_e - tss
        else:
            rel_s = tss - cand_e
            rel_e = tss - cand_s
        i_start = rel_s + window
        i_end = rel_e + window
        np.clip(i_start, 0, length, out=i_start)
        np.clip(i_end, 0, length, out=i_end)
        valid = i_end > i_start
        if not valid.any():
            continue
        np.add.at(diff, i_start[valid], 1)
        np.add.at(diff, i_end[valid], -1)
    arr = np.cumsum(diff[:-1])
    return arr, n_genes

# Statistical tests
def bh_correction(p_values):
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = np.array(p_values)[sorted_idx]
    adjusted = np.zeros(n)
    cummin = 1.0
    for i in range(n-1, -1, -1):
        rank = i + 1
        cummin = min(cummin, sorted_p[i] * n / rank)
        adjusted[sorted_idx[i]] = cummin
    return adjusted

def bootstrap_fold_diff_test(curves_h, curves_r, positions,
                             bin_size=50, n_bootstrap=1000, seed=42):
    n_pos = len(positions)
    n_bins = n_pos // bin_size
    bin_edges = np.linspace(0, n_pos, n_bins + 1, dtype=int)
    n_h = curves_h.shape[0]; n_r = curves_r.shape[0]
    rng = np.random.default_rng(seed)
    results = []
    for i in range(n_bins):
        start, end = bin_edges[i], bin_edges[i+1]
        h_vals = curves_h[:, start:end].mean(axis=1)
        r_vals = curves_r[:, start:end].mean(axis=1)
        obs_diff = h_vals.mean() - r_vals.mean()
        boot_diffs = np.zeros(n_bootstrap)
        for b in range(n_bootstrap):
            h_boot = rng.choice(h_vals, n_h, replace=True).mean()
            r_boot = rng.choice(r_vals, n_r, replace=True).mean()
            boot_diffs[b] = h_boot - r_boot
        ci_low = float(np.percentile(boot_diffs, 2.5))
        ci_high = float(np.percentile(boot_diffs, 97.5))
        p_val = 2 * min((boot_diffs <= 0).mean(), (boot_diffs >= 0).mean())
        p_val = float(min(p_val, 1.0))
        results.append({
            'bin_start': positions[start],
            'bin_end': positions[end-1],
            'bin_center': (positions[start] + positions[end-1]) / 2,
            'human_mean': float(h_vals.mean()),
            'rice_mean': float(r_vals.mean()),
            'diff': float(obs_diff),
            'ci_low': ci_low, 'ci_high': ci_high,
            'p_value': p_val,
        })
    df = pd.DataFrame(results)
    df['p_adj'] = bh_correction(df['p_value'].values)
    df['significant'] = df['p_adj'] < 0.05
    return df

def bootstrap_confidence_interval(curves, n_bootstrap=1000, ci=95, seed=42):
    n_samples, n_positions = curves.shape
    boot_means = np.zeros((n_bootstrap, n_positions))
    rng = np.random.default_rng(seed)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n_samples, n_samples)
        boot_means[i] = curves[idx].mean(axis=0)
    lower = np.percentile(boot_means, (100 - ci) / 2, axis=0)
    upper = np.percentile(boot_means, 100 - (100 - ci) / 2, axis=0)
    return lower, upper

# Per-species processing
def process_species(species, genome_dir, peak_dir, outdir, window,
                    smooth_sigma, min_qvalue=None):
    genome_dir = Path(genome_dir); peak_dir = Path(peak_dir); outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    label = species.capitalize()
    chrom_sizes = read_chrom_sizes(genome_dir / f'{species}.chrom.sizes')
    genes = read_genes_bed(genome_dir / f'{species}_genes.bed')
    if genes.empty or not chrom_sizes:
        raise SystemExit(f'[{label}] 缺少 genes.bed 或 chrom.sizes')
    genome_total = sum(chrom_sizes.values())
    print(f'[{label}] genes={len(genes)}  genome={genome_total:,} bp')
    if species == 'human':
        display_map = HUMAN_DISPLAY; color_map = COLORS_HUMAN; groups = HUMAN_GROUPS
    else:
        display_map = RICE_DISPLAY; color_map = COLORS_RICE; groups = RICE_GROUPS
    peak_files = []
    for p in sorted(peak_dir.iterdir()):
        if not p.is_file():
            continue
        if not p.name.endswith(('.bed', '.bed.gz', '.narrowPeak', '.narrowPeak.gz')):
            continue
        sample_name = strip_suffix(p.name)
        if 'G4' in sample_name.upper():
            print(f'  [skip G4] {sample_name}')
            continue
        peak_files.append(p)
    if not peak_files:
        raise SystemExit(f'[{label}] 没有可用的非 G4 peak 文件')
    positions = np.arange(-window, window + 1)
    results = {}
    for pf in peak_files:
        sample = strip_suffix(pf.name)
        print(f'  [{sample}]', end='  ', flush=True)
        peaks = read_peak_bed(pf, min_qvalue=min_qvalue)
        if peaks.empty:
            print(f'空文件（过滤后），跳过（min_qvalue={min_qvalue}）')
            continue
        peaks = peaks[peaks['chrom'].isin(chrom_sizes.keys())]
        if peaks.empty:
            print('染色体名不匹配，跳过')
            continue
        lookup = build_peak_lookup(peaks)
        arr, n_genes = accumulate_tss(genes, lookup, window)
        peak_bp = merged_length(peaks)
        peak_density = peak_bp / genome_total if genome_total > 0 else 0.0
        density = arr / n_genes if n_genes > 0 else np.zeros(len(positions))
        fold = density / peak_density if peak_density > 0 else np.zeros(len(positions))
        fold_smooth = gaussian_filter1d(fold, sigma=smooth_sigma) if smooth_sigma > 0 else fold
        pd.DataFrame({
            'position': positions, 'density_per_gene': density, 'fold_vs_genome': fold,
        }).to_csv(outdir / f'{sample}_TSS.tsv', sep='\t', index=False)
        results[sample] = {
            'species': species, 'display': display_map.get(sample, sample),
            'color': color_map.get(sample, '#333333'),
            'n_peaks': len(peaks), 'n_genes': n_genes,
            'peak_bp': peak_bp, 'peak_density': peak_density,
            'fold_smooth': fold_smooth, 'positions': positions,
        }
        print(f'peaks={len(peaks):>7d}  peak_density={peak_density:.6f}  '
              f'max_fold={fold_smooth.max():.2f}')
    return results, groups

# Per-species grid
def plot_species_grid(label, results, outdir, window):
    samples = list(results.keys()); n = len(samples)
    if n == 0:
        return
    ncols = min(2, n); nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 5.0 * nrows),
                             squeeze=False)
    for i, s in enumerate(samples):
        ax = axes[i // ncols][i % ncols]; res = results[s]
        ax.fill_between(res['positions'], 1, res['fold_smooth'],
                        color=res['color'], alpha=0.15, zorder=2)
        ax.plot(res['positions'], res['fold_smooth'],
                color=res['color'], lw=1.8, zorder=3)
        ax.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
        ax.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
        ax.set_xlabel('Distance from TSS (bp)', fontsize=11)
        ax.set_ylabel('Fold enrichment vs genome', fontsize=11)
        ax.set_xticks(np.arange(-window, window + 1, window // 4))
        ax.set_xlim(-window, window)
        ax.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(res['display'], fontsize=12, fontweight='bold', loc='left', pad=8)
        text = (f"peaks = {res['n_peaks']:,}\nmax fold = {res['fold_smooth'].max():.2f}")
        ax.text(0.97, 0.97, text, transform=ax.transAxes, ha='right', va='top',
                fontsize=9.5,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          edgecolor='#CCCCCC', linewidth=0.7), zorder=10)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis('off')
    plt.tight_layout()
    fig.savefig(outdir / f'grid_TSS_{label.lower()}.png', dpi=300, bbox_inches='tight')
    fig.savefig(outdir / f'grid_TSS_{label.lower()}.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {outdir}/grid_TSS_{label.lower()}.png')

# Per-species group mean
def plot_group_mean(label, results, groups, outdir, window):
    if not groups:
        return
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for gname, gsamples, gcolor in groups:
        valid = [s for s in gsamples if s in results]
        if not valid:
            continue
        pos = results[valid[0]]['positions']
        curves = np.vstack([results[s]['fold_smooth'] for s in valid])
        n = curves.shape[0]; m = curves.mean(axis=0)
        if n > 1:
            sem = curves.std(axis=0, ddof=1) / np.sqrt(n)
            ax.fill_between(pos, np.maximum(m - sem, 1e-3), m + sem,
                            color=gcolor, alpha=0.20, zorder=2)
        ax.plot(pos, m, color=gcolor, lw=2.2, zorder=3, label=f'{gname} (n={n})')
    ax.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=11)
    ax.set_ylabel('Fold enrichment vs genome', fontsize=11)
    ax.set_xticks(np.arange(-window, window + 1, window // 4))
    ax.set_xlim(-window, window)
    ax.legend(loc='upper right', frameon=False, fontsize=10)
    ax.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(f'{label}: group-mean TSS enrichment (mean ± SEM)',
                 fontsize=12, fontweight='bold', loc='left', pad=10)
    plt.tight_layout()
    fig.savefig(outdir / f'group_mean_TSS_{label.lower()}.png',
                dpi=300, bbox_inches='tight')
    fig.savefig(outdir / f'group_mean_TSS_{label.lower()}.pdf',
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {outdir}/group_mean_TSS_{label.lower()}.png')

# Combined plot 1: facets (main figure, no CI, no stats box)
def plot_combined_facet(results_h, results_r, outdir, window):
    """
    单样本曲线图：只画曲线 + 曲线下方填充装饰。
    不画 CI（CI 放均值图），不画统计框。
    """
    fig = plt.figure(figsize=(10, 8.5))
    gs = GridSpec(2, 1, hspace=0.30,
                  left=0.10, right=0.96, top=0.94, bottom=0.09)
    ax_h = fig.add_subplot(gs[0])
    ax_r = fig.add_subplot(gs[1], sharex=ax_h)

    # Top: human
    for s, res in results_h.items():
        ax_h.fill_between(res['positions'], 1, res['fold_smooth'],
                          color=res['color'], alpha=0.15, zorder=2)
        ax_h.plot(res['positions'], res['fold_smooth'],
                  color=res['color'], lw=2.0, zorder=3,
                  label=res['display'])

    ax_h.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax_h.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax_h.set_ylabel('Fold enrichment\nvs genome', fontsize=11)
    ax_h.set_ylim(bottom=0)
    ax_h.set_xticks(np.arange(-window, window + 1, window // 4))
    ax_h.tick_params(labelbottom=False)
    ax_h.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax_h.set_axisbelow(True)
    ax_h.legend(loc='upper right', frameon=False, fontsize=9.5)
    ax_h.text(-0.08, 1.02, 'Human CUT&Tag',
              transform=ax_h.transAxes,
              fontsize=12, fontweight='bold', va='bottom', ha='left')

    # Bottom: rice
    for s, res in results_r.items():
        ax_r.fill_between(res['positions'], 1, res['fold_smooth'],
                          color=res['color'], alpha=0.15, zorder=2)
        ax_r.plot(res['positions'], res['fold_smooth'],
                  color=res['color'], lw=2.0, zorder=3,
                  label=res['display'])

    ax_r.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax_r.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax_r.set_xlabel('Distance from TSS (bp)', fontsize=11.5)
    ax_r.set_ylabel('Fold enrichment\nvs genome', fontsize=11)
    ax_r.set_ylim(bottom=0)
    ax_r.set_xticks(np.arange(-window, window + 1, window // 4))
    ax_r.set_xlim(-window, window)
    ax_r.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax_r.set_axisbelow(True)
    handles, labels = ax_r.get_legend_handles_labels()
    if len(handles) > 1:
        order = ([i for i, l in enumerate(labels) if 'seedling' in l.lower()]
                 + [i for i, l in enumerate(labels) if 'seedling' not in l.lower()])
        handles = [handles[i] for i in order]
        labels = [labels[i] for i in order]
    ax_r.legend(handles, labels, loc='upper right', frameon=False, fontsize=9.5)
    ax_r.text(-0.08, 1.02, 'Rice iM-IP',
              transform=ax_r.transAxes,
              fontsize=12, fontweight='bold', va='bottom', ha='left')

    fn_png = Path(outdir) / 'combined_TSS_facet.png'
    fn_pdf = Path(outdir) / 'combined_TSS_facet.pdf'
    fig.savefig(fn_png, dpi=300, bbox_inches='tight')
    fig.savefig(fn_pdf, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {fn_png}')

# Combined plot 2: log y-axis
def plot_combined_overlay_log(results_h, results_r, outdir, window):
    fig, ax = plt.subplots(figsize=(10, 6.5))
    for s, res in results_h.items():
        ax.plot(res['positions'], res['fold_smooth'], color=res['color'],
                lw=1.8, alpha=0.9, label=res['display'], zorder=3)
    for s, res in results_r.items():
        ax.plot(res['positions'], res['fold_smooth'], color=res['color'],
                lw=1.8, alpha=0.9, label=res['display'], zorder=3)
    ax.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=11.5)
    ax.set_ylabel('Fold enrichment vs genome (log scale)', fontsize=11.5)
    ax.set_xticks(np.arange(-window, window + 1, window // 4))
    ax.set_xlim(-window, window)
    ax.set_yscale('log'); ax.set_ylim(bottom=0.5)
    ax.legend(loc='upper right', frameon=False, fontsize=9.5, ncol=2)
    ax.grid(alpha=0.18, ls='--', lw=0.6, which='both', zorder=0)
    ax.set_axisbelow(True)
    ax.set_title('TSS enrichment: Human CUT&Tag vs Rice iM-IP (log y)',
                 fontsize=12, fontweight='bold', loc='left', pad=10)
    plt.tight_layout()
    fig.savefig(Path(outdir) / 'combined_TSS_overlay_log.png', dpi=300, bbox_inches='tight')
    fig.savefig(Path(outdir) / 'combined_TSS_overlay_log.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {outdir}/combined_TSS_overlay_log.png')

# Combined plot 3: normalized
def plot_combined_overlay_normalized(results_h, results_r, outdir, window):
    fig, ax = plt.subplots(figsize=(10, 6.5))
    for s, res in results_h.items():
        norm = res['fold_smooth'] / res['fold_smooth'].max()
        ax.plot(res['positions'], norm, color=res['color'],
                lw=1.8, alpha=0.9, label=res['display'], zorder=3)
    for s, res in results_r.items():
        norm = res['fold_smooth'] / res['fold_smooth'].max()
        ax.plot(res['positions'], norm, color=res['color'],
                lw=1.8, alpha=0.9, label=res['display'], zorder=3)
    ax.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax.axhline(0, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=11.5)
    ax.set_ylabel('Normalized enrichment\n(each curve / its max)', fontsize=11)
    ax.set_xticks(np.arange(-window, window + 1, window // 4))
    ax.set_xlim(-window, window)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='upper right', frameon=False, fontsize=9.5, ncol=2)
    ax.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title('TSS enrichment: normalized (shape comparison)',
                 fontsize=12, fontweight='bold', loc='left', pad=10)
    plt.tight_layout()
    fig.savefig(Path(outdir) / 'combined_TSS_overlay_normalized.png',
                dpi=300, bbox_inches='tight')
    fig.savefig(Path(outdir) / 'combined_TSS_overlay_normalized.pdf',
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {outdir}/combined_TSS_overlay_normalized.png')

# Combined plot 4: group-mean facets (supplementary A/B + stats box)
def plot_combined_species_mean_facet(results_h, results_r, outdir, window,
                                     bin_results=None,
                                     ks_stat=None, ks_p=None):
    fig = plt.figure(figsize=(9, 8.5))
    gs = GridSpec(2, 1, hspace=0.30,
                  left=0.10, right=0.96, top=0.94, bottom=0.09)
    ax_h = fig.add_subplot(gs[0])
    ax_r = fig.add_subplot(gs[1], sharex=ax_h)

    for gname, gsamples, gcolor in HUMAN_GROUPS:
        valid = [s for s in gsamples if s in results_h]
        if not valid:
            continue
        pos = results_h[valid[0]]['positions']
        curves = np.vstack([results_h[s]['fold_smooth'] for s in valid])
        n = curves.shape[0]; m = curves.mean(axis=0)
        if n > 1:
            sem = curves.std(axis=0, ddof=1) / np.sqrt(n)
            ax_h.fill_between(pos, np.maximum(m - sem, 0), m + sem,
                              color=gcolor, alpha=0.20, zorder=2)
        ax_h.plot(pos, m, color=gcolor, lw=2.6, zorder=3,
                  label=f'{gname} (n={n})')

    ax_h.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax_h.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax_h.set_ylabel('Fold enrichment\nvs genome', fontsize=11)
    ax_h.set_ylim(bottom=0)
    ax_h.set_xticks(np.arange(-window, window + 1, window // 4))
    ax_h.tick_params(labelbottom=False)
    ax_h.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax_h.set_axisbelow(True)
    ax_h.legend(loc='upper right', frameon=False, fontsize=10)
    ax_h.text(-0.08, 1.02, 'A', transform=ax_h.transAxes,
              fontsize=17, fontweight='bold', va='bottom', ha='left')
    ax_h.text(0.01, 1.02, 'Human CUT&Tag (mean ± SEM)',
              transform=ax_h.transAxes,
              fontsize=12, fontweight='bold', va='bottom', ha='left')

    for gname, gsamples, gcolor in RICE_GROUPS:
        valid = [s for s in gsamples if s in results_r]
        if not valid:
            continue
        pos = results_r[valid[0]]['positions']
        curves = np.vstack([results_r[s]['fold_smooth'] for s in valid])
        n = curves.shape[0]; m = curves.mean(axis=0)
        if n > 1:
            sem = curves.std(axis=0, ddof=1) / np.sqrt(n)
            ax_r.fill_between(pos, np.maximum(m - sem, 0), m + sem,
                              color=gcolor, alpha=0.20, zorder=2)
        ax_r.plot(pos, m, color=gcolor, lw=2.6, zorder=3,
                  label=f'{gname} (n={n})')

    ax_r.axvline(0, color='#333333', ls='--', lw=0.9, alpha=0.7, zorder=1)
    ax_r.axhline(1, color='#999999', ls=':', lw=0.8, alpha=0.6, zorder=1)
    ax_r.set_xlabel('Distance from TSS (bp)', fontsize=11.5)
    ax_r.set_ylabel('Fold enrichment\nvs genome', fontsize=11)
    ax_r.set_ylim(bottom=0)
    ax_r.set_xticks(np.arange(-window, window + 1, window // 4))
    ax_r.set_xlim(-window, window)
    ax_r.grid(alpha=0.18, ls='--', lw=0.6, zorder=0)
    ax_r.set_axisbelow(True)
    ax_r.legend(loc='upper right', frameon=False, fontsize=10)
    ax_r.text(-0.08, 1.02, 'B', transform=ax_r.transAxes,
              fontsize=17, fontweight='bold', va='bottom', ha='left')
    ax_r.text(0.01, 1.02, 'Rice iM-IP (mean ± SEM)',
              transform=ax_r.transAxes,
              fontsize=12, fontweight='bold', va='bottom', ha='left')

    # Stats box
    if bin_results is not None:
        n_sig = int(bin_results['significant'].sum())
        n_total = len(bin_results)

        stat_lines = [
            f"Bootstrap: {n_sig}/{n_total} bins significant (FDR < 0.05)",
        ]
        if ks_stat is not None and ks_p is not None:
            p_str = "P < 0.001" if ks_p < 0.001 else f"P = {ks_p:.3e}"
            stat_lines.append(
                f"Kolmogorov-Smirnov: D = {ks_stat:.3f}, {p_str}"
            )

        stat_text = '\n'.join(stat_lines)

        ax_h.text(0.98, 0.62, stat_text,
                  transform=ax_h.transAxes,
                  ha='right', va='top',
                  fontsize=8.5,
                  bbox=dict(boxstyle='round,pad=0.4',
                            facecolor='#FAFAFA',
                            edgecolor='#CCCCCC',
                            linewidth=0.7),
                  zorder=10)

    fig.savefig(Path(outdir) / 'combined_TSS_species_mean_facet.png',
                dpi=300, bbox_inches='tight')
    fig.savefig(Path(outdir) / 'combined_TSS_species_mean_facet.pdf',
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved: {outdir}/combined_TSS_species_mean_facet.png')

# main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--human_genome_dir', required=True)
    ap.add_argument('--human_peak_dir',   required=True)
    ap.add_argument('--rice_genome_dir',  required=True)
    ap.add_argument('--rice_peak_dir',    required=True)
    ap.add_argument('--outdir',           required=True)
    ap.add_argument('--window',       type=int,   default=1000)
    ap.add_argument('--smooth_sigma', type=float, default=30)
    ap.add_argument('--rice_min_qvalue', type=float, default=10.0)
    ap.add_argument('--bin_size', type=int, default=50)
    ap.add_argument('--n_bootstrap', type=int, default=1000)
    args = ap.parse_args()

    set_nature_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print('=' * 60)
    print('Processing HUMAN CUT&Tag (hg38)...')
    human_results, human_groups = process_species(
        'human', args.human_genome_dir, args.human_peak_dir,
        outdir / 'human', args.window, args.smooth_sigma, min_qvalue=None)

    print('=' * 60)
    print(f'Processing RICE (min_qvalue={args.rice_min_qvalue})...')
    rice_results, rice_groups = process_species(
        'rice', args.rice_genome_dir, args.rice_peak_dir,
        outdir / 'rice', args.window, args.smooth_sigma,
        min_qvalue=args.rice_min_qvalue)

    print('=' * 60)
    print('Performing statistical tests...')
    curves_h = np.vstack([res['fold_smooth'] for res in human_results.values()])
    curves_r = np.vstack([res['fold_smooth'] for res in rice_results.values()])
    positions = list(human_results.values())[0]['positions']

    bin_df = bootstrap_fold_diff_test(curves_h, curves_r, positions,
                                       bin_size=args.bin_size,
                                       n_bootstrap=args.n_bootstrap)
    bin_df.to_csv(outdir / 'bootstrap_binned_results.tsv', sep='\t', index=False)
    n_sig = bin_df['significant'].sum()
    print(f'  Bootstrap: {n_sig}/{len(bin_df)} bins significant (FDR < 0.05)')

    # CI is still computed for possible later use, but the main figure does not use it
    lower_h, upper_h = bootstrap_confidence_interval(curves_h, n_bootstrap=args.n_bootstrap)
    lower_r, upper_r = bootstrap_confidence_interval(curves_r, n_bootstrap=args.n_bootstrap)

    mean_h = curves_h.mean(axis=0); mean_r = curves_r.mean(axis=0)
    ks_stat, ks_p = ks_2samp(mean_h, mean_r)
    print(f'  Kolmogorov-Smirnov: D = {ks_stat:.3f}, p = {ks_p:.3e}')

    with open(outdir / 'statistical_summary.txt', 'w') as f:
        f.write(f'Bootstrap significant bins: {n_sig}/{len(bin_df)}\n')
        f.write(f'Kolmogorov-Smirnov: D = {ks_stat:.4f}, p = {ks_p:.3e}\n')

    print('=' * 60)
    print('Generating per-species plots...')
    plot_species_grid('Human', human_results, outdir / 'human', args.window)
    plot_group_mean('Human', human_results, human_groups, outdir / 'human', args.window)
    plot_species_grid('Rice', rice_results, outdir / 'rice', args.window)
    plot_group_mean('Rice', rice_results, rice_groups, outdir / 'rice', args.window)

    print('=' * 60)
    print('Generating combined plots...')
    # Main figure: no CI, no stats box
    plot_combined_facet(human_results, rice_results, outdir, args.window)
    # Supplementary figure: mean + SEM + stats box
    plot_combined_species_mean_facet(human_results, rice_results, outdir, args.window,
                                     bin_results=bin_df,
                                     ks_stat=ks_stat, ks_p=ks_p)
    plot_combined_overlay_log(human_results, rice_results, outdir, args.window)
    plot_combined_overlay_normalized(human_results, rice_results, outdir, args.window)

    print(f'\nAll done! Results in {outdir}')

if __name__ == '__main__':
    main()
