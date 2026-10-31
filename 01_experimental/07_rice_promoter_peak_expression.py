#!/usr/bin/env python3
"""
011_rice_promoter_peak_expression.py  (v2 - Nature color scheme + editable PDF)

Rice: expression distribution of genes with a seedling IP peak in the promoter vs all genes

- Statistical test: full data, one-sided Mann-Whitney U
- Plotting: after IQR outlier removal
- Three plots per comparison: violin / bar / box
- Each comparison split by: raw TPM / log10(TPM+1)
- Editable PDF (pdf.fonttype = 42)

Colors:
  iM IP-seq seedling peak near TSS : dark green #238B45
  All Genes                        : neutral gray #BDBDBD
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

np.random.seed(42)


# Global Nature style + editable PDF
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
COLOR_PEAK = '#238B45'   # dark green (Rice iM IP-seq)
COLOR_ALL  = '#BDBDBD'   # neutral gray (All Genes)


# Configuration
PROMOTER = "input/promoters_TSSpm1kb.rice.bed"
PEAK_SEEDLING = "input/IP_vs_Input_intersect.narrowPeak"
TPM_FILE = "input/NP_mock_TPM_Os.csv"
OUT_TABLE = "rice_gene_expression_peak_status.tsv"


# Find genes with a peak in the promoter region
def get_genes_with_peak(promoter, peak, label):
    print(f"--- Find genes with a peak in the promoter region: {label} ---")
    cmd = (f"grep -v '^#' {promoter} | "
           f"bedtools intersect -a - -b {peak} -wa -u")
    result = subprocess.run(cmd, shell=True,
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bedtools failed ({label}):\n{result.stderr}")
    genes = set()
    for line in result.stdout.strip().split('\n'):
        if line.strip() == '':
            continue
        parts = line.split('\t')
        if len(parts) >= 4:
            genes.add(parts[3])
    print(f"    Genes with a peak: {len(genes)}")
    return genes


print("=== 1. Intersect promoters with peaks ===")
genes_seedling = get_genes_with_peak(PROMOTER, PEAK_SEEDLING,
                                     "Rice IP-seq (seedling)")

# Read TPM
print("\n=== 2. Read TPM ===")
tpm = pd.read_csv(TPM_FILE)
tpm_cols = [c for c in tpm.columns if c.startswith('NP_mock')]
print(f"TPM replicate columns: {tpm_cols}")

tpm['TPM_mean'] = tpm[tpm_cols].mean(axis=1)
tpm_exp = tpm[(tpm[tpm_cols] > 0).any(axis=1)].copy()
print(f"Expressed genes: {len(tpm_exp)}")

tpm_exp['has_seedling'] = tpm_exp['gene_id'].isin(genes_seedling)

n_seedling = int(tpm_exp['has_seedling'].sum())
n_all = len(tpm_exp)
print(f"Genes with a seedling peak and expression: {n_seedling}")
print(f"All expressed genes: {n_all}")

# Extract expression values
seedling_expr = tpm_exp[tpm_exp['has_seedling']]['TPM_mean'].values
all_expr = tpm_exp['TPM_mean'].values

# Statistical test (full data)
def p_to_star(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'


_, p_seed = stats.mannwhitneyu(seedling_expr, all_expr,
                               alternative='greater')
star_seed = p_to_star(p_seed)
print(f"\nseedling vs All: p = {p_seed:.3e} ({star_seed}), "
      f"median seedling={np.median(seedling_expr):.3f}, "
      f"All={np.median(all_expr):.3f}")

# IQR outlier removal
def remove_outliers(x):
    if len(x) == 0:
        return x
    q25, q75 = np.percentile(x, [25, 75])
    iqr = q75 - q25
    lower = q25 - 1.5 * iqr
    upper = q75 + 1.5 * iqr
    return x[(x >= lower) & (x <= upper)]


seedling_clean = remove_outliers(seedling_expr)
all_clean = remove_outliers(all_expr)
print(f"After outlier removal: seedling={len(seedling_clean)}, All={len(all_clean)}")


# Shared plotting helpers
def _prepare(data_list, log_transform):
    if log_transform:
        return [np.log10(x + 1) for x in data_list], 'log10(TPM + 1)'
    return data_list, 'Gene expression (TPM)'


def _sig_bar(ax, x1, x2, y, star, h):
    ax.plot([x1, x1, x2, x2],
            [y, y + h, y + h, y],
            color='black', lw=1.0, zorder=10)
    ax.text((x1 + x2) / 2, y + h * 1.1, star,
            ha='center', va='bottom',
            fontsize=13, fontweight='bold', zorder=10)


def _draw_jitter(ax, data_list, colors, max_points=1500,
                 point_size=1.4, alpha=0.18):
    for i, data in enumerate(data_list):
        x_pos = i + 1
        if len(data) > max_points:
            data_show = np.random.choice(data, max_points, replace=False)
        else:
            data_show = data
        jitter = np.random.normal(0, 0.06, size=len(data_show))
        ax.scatter(x_pos + jitter, data_show,
                   s=point_size, alpha=alpha, color=colors[i],
                   edgecolor='none', rasterized=True, zorder=2)


# Figure 1: violin + jittered scatter points
def draw_violin(data_list, labels, colors, stars,
                out_prefix, title,
                log_transform=False, max_points=1500):
    data_plot_list, ylabel = _prepare(data_list, log_transform)

    fig, ax = plt.subplots(figsize=(4.8, 5.0))
    positions = list(range(1, len(data_plot_list) + 1))

    vp = ax.violinplot(data_plot_list, positions=positions,
                       widths=0.65, showmeans=False,
                       showmedians=False, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_alpha(0.55)
        body.set_edgecolor('white')
        body.set_linewidth(0.8)

    _draw_jitter(ax, data_plot_list, colors, max_points)

    y_max = max([np.max(d) for d in data_plot_list if len(d) > 0])
    y_base = y_max * 1.03
    h = y_max * 0.025
    _sig_bar(ax, 1, 2, y_base, stars[0], h)
    ax.set_ylim(top=y_max * 1.15)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=10.5)
    ax.set_ylabel(ylabel, fontsize=11.5)
    ax.set_ylim(bottom=0)
    ax.grid(axis='y', ls='--', lw=0.6, alpha=0.25)
    ax.set_axisbelow(True)

    ax.set_title(title, fontsize=12, fontweight='bold', loc='left', pad=8)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Violin: {out_prefix}.pdf/.png")


# Figure 2: bar chart (mean ± SEM)
def draw_bar(data_list, labels, colors, stars,
             out_prefix, title,
             log_transform=False, bar_width=0.45):
    data_plot_list, _ = _prepare(data_list, log_transform)
    if log_transform:
        ylabel = 'Mean log10(TPM + 1)'
    else:
        ylabel = 'Mean gene expression (TPM)'

    means = [np.mean(x) for x in data_plot_list]
    sems  = [np.std(x, ddof=1) / np.sqrt(len(x)) for x in data_plot_list]

    fig, ax = plt.subplots(figsize=(4.0, 5.0))
    x_pos = np.arange(len(data_plot_list))

    ax.bar(x_pos, means, yerr=sems, width=bar_width,
           color=colors, alpha=0.9,
           edgecolor='black', linewidth=0.8,
           error_kw=dict(ecolor='black', elinewidth=1.2,
                         capsize=5, capthick=1.2),
           zorder=3)

    y_top = max(m + s for m, s in zip(means, sems))
    y_base = y_top * 1.03
    h = y_top * 0.02
    _sig_bar(ax, 0, 1, y_base, stars[0], h)
    ax.set_ylim(top=y_top * 1.15)

    for x, m in zip(x_pos, means):
        ax.text(x, m * 0.5, f'{m:.2f}',
                ha='center', va='center',
                fontsize=10, fontweight='bold', color='white',
                zorder=4)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=10.5)
    ax.set_ylabel(ylabel, fontsize=11.5)
    ax.set_ylim(bottom=0)
    ax.grid(axis='y', ls='--', lw=0.6, alpha=0.25)
    ax.set_axisbelow(True)

    ax.set_title(title, fontsize=12, fontweight='bold', loc='left', pad=8)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Bar: {out_prefix}.pdf/.png")


# Figure 3: box plot + jittered scatter points
def draw_box(data_list, labels, colors, stars,
             out_prefix, title,
             log_transform=False, max_points=1500):
    data_plot_list, ylabel = _prepare(data_list, log_transform)

    fig, ax = plt.subplots(figsize=(4.8, 5.0))
    positions = list(range(1, len(data_plot_list) + 1))

    bp = ax.boxplot(data_plot_list, positions=positions,
                    widths=0.5, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(linewidth=1.2),
                    medianprops=dict(color='black', linewidth=1.6),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2))

    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.75)
        patch.set_edgecolor('black')
        patch.set_linewidth(1.0)

    _draw_jitter(ax, data_plot_list, colors, max_points)

    y_max = max([np.max(d) for d in data_plot_list if len(d) > 0])
    y_base = y_max * 1.03
    h = y_max * 0.025
    _sig_bar(ax, 1, 2, y_base, stars[0], h)
    ax.set_ylim(top=y_max * 1.15)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=10.5)
    ax.set_ylabel(ylabel, fontsize=11.5)
    ax.set_ylim(bottom=0)
    ax.grid(axis='y', ls='--', lw=0.6, alpha=0.25)
    ax.set_axisbelow(True)

    ax.set_title(title, fontsize=12, fontweight='bold', loc='left', pad=8)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Box: {out_prefix}.pdf/.png")


# Main workflow
set_nature_style()

labels = ['Rice iM-IP\n(seedling)', 'All Genes']
colors = [COLOR_PEAK, COLOR_ALL]
stars  = [star_seed, 'ns']

print("\n--- IQR-filtered + raw TPM ---")
draw_violin([seedling_clean, all_clean], labels, colors, stars,
            out_prefix="rice_iqr_raw_violin",
            title="Rice: IQR-filtered, raw TPM",
            log_transform=False)

draw_bar([seedling_clean, all_clean], labels, colors, stars,
         out_prefix="rice_iqr_raw_bar",
         title="Rice: IQR-filtered, mean TPM",
         log_transform=False)

draw_box([seedling_clean, all_clean], labels, colors, stars,
         out_prefix="rice_iqr_raw_box",
         title="Rice: IQR-filtered, raw TPM (box)",
         log_transform=False)

print("\n--- IQR-filtered + log10(TPM+1) ---")
draw_violin([seedling_clean, all_clean], labels, colors, stars,
            out_prefix="rice_iqr_log_violin",
            title="Rice: IQR-filtered, log10(TPM+1)",
            log_transform=True)

draw_bar([seedling_clean, all_clean], labels, colors, stars,
         out_prefix="rice_iqr_log_bar",
         title="Rice: IQR-filtered, mean log10(TPM+1)",
         log_transform=True)

draw_box([seedling_clean, all_clean], labels, colors, stars,
         out_prefix="rice_iqr_log_box",
         title="Rice: IQR-filtered, log10(TPM+1) (box)",
         log_transform=True)

# Save table
tpm_exp[['gene_id', 'TPM_mean', 'has_seedling']].to_csv(
    OUT_TABLE, sep='\t', index=False)
print(f"\nSaved as {OUT_TABLE}")
print("All done.")
