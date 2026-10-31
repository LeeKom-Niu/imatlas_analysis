#!/usr/bin/env python3
"""
022_plot_hek293t_promoter_expression.py  (v2 - Nature 配色 + 可编辑 PDF)

HEK293T：CUT&Tag iM 启动子关联基因 vs 全部基因的表达分布

- 统计检验：完整数据，Mann-Whitney U 单侧
- 绘图：IQR 去离群后
- 每个组合画三张图：小提琴 / 柱状图 / 箱线图
- 每个组合又分：原始 TPM / log10(TPM+1)
- 可编辑 PDF（pdf.fonttype = 42）

配色：
  iM CUT&Tag near TSS : 深蓝 #08519C
  All Genes           : 中性灰 #BDBDBD
"""

import re
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
COLOR_PEAK = '#08519C'   # dark blue (Human iM CUT&Tag)
COLOR_ALL  = '#BDBDBD'   # neutral gray (All Genes)


# 1. Read expression data
gene_tpm = pd.read_csv("input/HEK293T_gene_TPM.tsv", sep='\t')

gene_tpm_exp = gene_tpm[
    (gene_tpm['TPM_rep1'] > 0) |
    (gene_tpm['TPM_rep2'] > 0) |
    (gene_tpm['TPM_rep3'] > 0)
].copy()


def strip_version(gene_id):
    return re.sub(r'\.\d+$', '', str(gene_id))


gene_tpm_exp['gene_id_clean'] = gene_tpm_exp['gene_id'].apply(strip_version)

# 2. Read the CUT&Tag overlap file
cut_tag_overlap = pd.read_csv("CUTTag_iM_promoter_overlap.bed",
                              sep='\t', header=None)
cut_tag_genes = set(
    cut_tag_overlap.iloc[:, 9].apply(strip_version).unique()
)
print(f"CUT&Tag 启动子关联基因数: {len(cut_tag_genes)}")

# 3. Extract expression values
cut_tag_expr = gene_tpm_exp[
    gene_tpm_exp['gene_id_clean'].isin(cut_tag_genes)
]['TPM_mean'].values
all_expr = gene_tpm_exp['TPM_mean'].values

print(f"CUT&Tag 匹配到表达的基因数: {len(cut_tag_expr)}")
print(f"全部表达基因数: {len(all_expr)}")

# 4. Statistical test (full data)
def p_to_star(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'


_, p_cut = stats.mannwhitneyu(cut_tag_expr, all_expr, alternative='greater')
star_cut = p_to_star(p_cut)
print(f"CUT&Tag vs All: p = {p_cut:.3e} ({star_cut})")

# 5. IQR outlier removal
def remove_outliers(x):
    if len(x) == 0:
        return x
    q25, q75 = np.percentile(x, [25, 75])
    iqr = q75 - q25
    lower = q25 - 1.5 * iqr
    upper = q75 + 1.5 * iqr
    return x[(x >= lower) & (x <= upper)]


cut_tag_clean = remove_outliers(cut_tag_expr)
all_clean = remove_outliers(all_expr)
print(f"去离群后：CUT&Tag={len(cut_tag_clean)}, All={len(all_clean)}")


# Shared plotting helpers
def _prepare(data_list, log_transform):
    if log_transform:
        return [np.log10(x + 1) for x in data_list], 'log10(TPM + 1)'
    return data_list, 'Gene expression (TPM)'


def _sig_bar(ax, x1, x2, y, star, h):
    """显著性横线 + 星号"""
    ax.plot([x1, x1, x2, x2],
            [y, y + h, y + h, y],
            color='black', lw=1.0, zorder=10)
    ax.text((x1 + x2) / 2, y + h * 1.1, star,
            ha='center', va='bottom',
            fontsize=13, fontweight='bold', zorder=10)


def _draw_jitter(ax, data_list, colors, max_points=1500,
                 point_size=1.4, alpha=0.18):
    """抖散点"""
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

    # Violin
    vp = ax.violinplot(data_plot_list, positions=positions,
                       widths=0.65, showmeans=False,
                       showmedians=False, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_alpha(0.55)
        body.set_edgecolor('white')
        body.set_linewidth(0.8)

    # Jittered scatter points
    _draw_jitter(ax, data_plot_list, colors, max_points)

    # Significance bar
    y_max = max([np.max(d) for d in data_plot_list if len(d) > 0])
    y_base = y_max * 1.03
    h = y_max * 0.025
    _sig_bar(ax, 1, 2, y_base, stars[0], h)
    ax.set_ylim(top=y_max * 1.15)

    # Axes
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
    print(f"  小提琴: {out_prefix}.pdf/.png")


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

    # Bars
    ax.bar(x_pos, means, yerr=sems, width=bar_width,
           color=colors, alpha=0.9,
           edgecolor='black', linewidth=0.8,
           error_kw=dict(ecolor='black', elinewidth=1.2,
                         capsize=5, capthick=1.2),
           zorder=3)

    # Significance bar
    y_top = max(m + s for m, s in zip(means, sems))
    y_base = y_top * 1.03
    h = y_top * 0.02
    _sig_bar(ax, 0, 1, y_base, stars[0], h)
    ax.set_ylim(top=y_top * 1.15)

    # Value labels
    for x, m in zip(x_pos, means):
        ax.text(x, m * 0.5, f'{m:.2f}',
                ha='center', va='center',
                fontsize=10, fontweight='bold', color='white',
                zorder=4)

    # Axes
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
    print(f"  柱状图: {out_prefix}.pdf/.png")


# Figure 3: box plot + jittered scatter points
def draw_box(data_list, labels, colors, stars,
             out_prefix, title,
             log_transform=False, max_points=1500):
    data_plot_list, ylabel = _prepare(data_list, log_transform)

    fig, ax = plt.subplots(figsize=(4.8, 5.0))
    positions = list(range(1, len(data_plot_list) + 1))

    # Box plot
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

    # Jittered scatter points
    _draw_jitter(ax, data_plot_list, colors, max_points)

    # Significance bar
    y_max = max([np.max(d) for d in data_plot_list if len(d) > 0])
    y_base = y_max * 1.03
    h = y_max * 0.025
    _sig_bar(ax, 1, 2, y_base, stars[0], h)
    ax.set_ylim(top=y_max * 1.15)

    # Axes
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
    print(f"  箱线图: {out_prefix}.pdf/.png")


# Main workflow
set_nature_style()

labels = ['iM CUT&Tag\nnear TSS', 'All Genes']
colors = [COLOR_PEAK, COLOR_ALL]
stars  = [star_cut, 'ns']

print("\n--- 去离群 + 原始 TPM ---")
draw_violin([cut_tag_clean, all_clean], labels, colors, stars,
            out_prefix="hek293t_iqr_raw_violin",
            title="HEK293T: IQR-filtered, raw TPM",
            log_transform=False)

draw_bar([cut_tag_clean, all_clean], labels, colors, stars,
         out_prefix="hek293t_iqr_raw_bar",
         title="HEK293T: IQR-filtered, mean TPM",
         log_transform=False)

draw_box([cut_tag_clean, all_clean], labels, colors, stars,
         out_prefix="hek293t_iqr_raw_box",
         title="HEK293T: IQR-filtered, raw TPM (box)",
         log_transform=False)

print("\n--- 去离群 + log10(TPM+1) ---")
draw_violin([cut_tag_clean, all_clean], labels, colors, stars,
            out_prefix="hek293t_iqr_log_violin",
            title="HEK293T: IQR-filtered, log10(TPM+1)",
            log_transform=True)

draw_bar([cut_tag_clean, all_clean], labels, colors, stars,
         out_prefix="hek293t_iqr_log_bar",
         title="HEK293T: IQR-filtered, mean log10(TPM+1)",
         log_transform=True)

draw_box([cut_tag_clean, all_clean], labels, colors, stars,
         out_prefix="hek293t_iqr_log_box",
         title="HEK293T: IQR-filtered, log10(TPM+1) (box)",
         log_transform=True)

print("\n全部完成。")
