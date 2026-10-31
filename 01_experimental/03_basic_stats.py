#!/usr/bin/env python3
"""
011_basic_stats.py  (v3 - Nature 配色 + 可编辑 PDF + 顺序/颜色调整)

对所有 iM 数据（人类 + 水稻）进行基础统计，并将峰值数量绘制在同一柱形图中。
输出目录：011_basic_stats

人类数据：
  - HEK293T / WDLPS : CUT&Tag (iM)
  - GSE248445       : DIP-seq (iMab 抗体, 3 个单细胞 + 1 个 3-cell shared)
水稻数据：
  - IP_vs_IgG_intersect.narrowPeak    ← GSE234751 CK-7 iM-IP vs IgG（pH 7.0，未处理）
  - IP_vs_Input_intersect.narrowPeak  ← GSE184783 iM-IP vs Input（水稻幼苗）

Nature 配色：
  - Human CUT&Tag : 深蓝
  - Human DIP-seq : 浅蓝（同一色系）
  - Rice Input    : 深绿（seedling，重点）
  - Rice IgG      : 浅绿（同一色系）
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')


# Global Nature style + editable PDF
mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'axes.linewidth': 0.9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10.5,
    'xtick.major.width': 0.9,
    'ytick.major.width': 0.9,
    'xtick.major.size': 3.5,
    'ytick.major.size': 3.5,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'legend.fontsize': 10,
    'legend.frameon': False,
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    # Editable PDF (TrueType)
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'pdf.compression': 6,
})


# Color scheme
# Human: blue palette (CUT&Tag dark, DIP-seq light)
# Rice: green palette (Input dark, IgG light) — seedling is the focus
COLOR_HUMAN_CUTTAG = '#08519C'   # dark blue
COLOR_HUMAN_DIPSEQ = '#9ECAE1'   # light blue

# Change 1: swap the two rice colors (seedling = Input uses dark green)
COLOR_RICE_INPUT   = '#238B45'   # dark green (previously used for IgG)
COLOR_RICE_IGG     = '#A1D99B'   # light green (previously used for Input)

# Exact mapping by sample name
SAMPLE_COLORS = {
    # Human CUT&Tag
    'HEK293T (CUT&Tag)':  COLOR_HUMAN_CUTTAG,
    'WDLPS (CUT&Tag)':    COLOR_HUMAN_CUTTAG,
    # Human DIP-seq
    'HEK293T (DIP-seq)':  COLOR_HUMAN_DIPSEQ,
    'MCF7 (DIP-seq)':     COLOR_HUMAN_DIPSEQ,
    'U2OS (DIP-seq)':     COLOR_HUMAN_DIPSEQ,
    '3-cell (DIP-seq (shared))': COLOR_HUMAN_DIPSEQ,
    # Rice
    'Rice pH7.0 (iM-IP vs IgG)':      COLOR_RICE_IGG,
    'Rice seedling (iM-IP vs Input)': COLOR_RICE_INPUT,
}


# Configuration
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/19_modified_analysis/10_iMab/05_Experiment_analysis"
DATA_DIR = os.path.join(BASE_DIR, "experiment_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "011_basic_stats")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# Scan data files
def scan_files():
    files = []

    # Human
    human_dir = os.path.join(DATA_DIR, "human")
    for f in sorted(os.listdir(human_dir)):
        if not (f.endswith('.bed') or f.endswith('.narrowPeak') or f.endswith('.bed.gz')):
            continue
        file_path = os.path.join(human_dir, f)

        if 'HEK293T' in f and 'iM_highConf' in f:
            cell = 'HEK293T'; method = 'CUT&Tag';            order = 1
        elif 'WDLPS' in f and 'iM_highConf' in f:
            cell = 'WDLPS';   method = 'CUT&Tag';            order = 2
        elif 'GSE248445' in f:
            if 'MCF7_HEK293T_U2OS' in f:
                cell = '3-cell';   method = 'DIP-seq (shared)'; order = 6
            # Change 2: swap the order of HEK293T DIP and MCF7 DIP
            elif 'HEK293T' in f and 'U2OS' not in f and 'MCF7' not in f:
                cell = 'HEK293T';  method = 'DIP-seq';          order = 4   # previously 3
            elif 'MCF7' in f and '16C' not in f and 'HEK293T' not in f and 'U2OS' not in f:
                cell = 'MCF7';     method = 'DIP-seq';          order = 3   # previously 4
            elif 'U2OS' in f and 'HEK293T' not in f:
                cell = 'U2OS';     method = 'DIP-seq';          order = 5
            else:
                cell = 'unknown';  method = 'DIP-seq';          order = 7
        else:
            continue

        files.append({
            'file': file_path,
            'name': f,
            'display_name': f"{cell} ({method})",
            'species': 'human',
            'cell': cell,
            'method': method,
            'order': order,
            'is_gz': f.endswith('.gz'),
        })

    # Rice
    rice_dir = os.path.join(DATA_DIR, "rice")
    for f in sorted(os.listdir(rice_dir)):
        if not (f.endswith('.narrowPeak') or f.endswith('.narrowPeak.gz')):
            continue
        file_path = os.path.join(rice_dir, f)

        if 'IgG' in f:
            cell = 'Rice pH7.0'
            method = 'iM-IP vs IgG'
            order = 1
        elif 'Input' in f:
            cell = 'Rice seedling'
            method = 'iM-IP vs Input'
            order = 2
        else:
            cell = 'Rice'
            method = 'iM-IP-seq'
            order = 3

        files.append({
            'file': file_path,
            'name': f,
            'display_name': f"{cell} ({method})",
            'species': 'rice',
            'cell': cell,
            'method': method,
            'order': order,
            'is_gz': f.endswith('.gz'),
        })

    files.sort(key=lambda x: (x['species'] == 'rice', x['order']))
    human_files = [f for f in files if f['species'] == 'human']
    rice_files  = [f for f in files if f['species'] == 'rice']
    for i, f in enumerate(human_files):
        f['order'] = i
    for i, f in enumerate(rice_files):
        f['order'] = i
    files = human_files + rice_files
    return files


# Read BED files and compute statistics
def read_bed_stats(file_info):
    file_path = file_info['file']
    is_gz = file_info['is_gz']
    if is_gz:
        df = pd.read_csv(file_path, sep='\t', header=None, compression='gzip')
    else:
        df = pd.read_csv(file_path, sep='\t', header=None)

    n_cols   = len(df.columns)
    n_peaks  = len(df)
    widths   = df.iloc[:, 2] - df.iloc[:, 1]
    width_mean   = widths.mean()
    width_median = widths.median()
    width_std    = widths.std()
    width_sem    = width_std / np.sqrt(n_peaks) if n_peaks > 0 else 0

    signal = None
    signal_col = None
    signal_mean = signal_median = signal_std = signal_sem = None

    if n_cols >= 6:
        try:
            signal = df.iloc[:, 4].astype(float)
            signal_col = 'max_signal'
            if signal.median() == 0 or signal.isna().all():
                signal = df.iloc[:, 3].astype(float)
                signal_col = 'AUC'
        except Exception:
            pass
    if n_cols >= 10 and signal is None:
        try:
            signal = df.iloc[:, 6].astype(float)
            signal_col = 'signalValue'
        except Exception:
            pass

    if signal is not None and len(signal) > 0:
        signal_mean   = signal.mean()
        signal_median = signal.median()
        signal_std    = signal.std()
        signal_sem    = signal_std / np.sqrt(len(signal))
    else:
        signal_col = None

    return {
        'name': file_info['name'],
        'display_name': file_info['display_name'],
        'species': file_info['species'],
        'cell': file_info['cell'],
        'method': file_info['method'],
        'order': file_info['order'],
        'n_peaks': n_peaks,
        'width_mean': width_mean,
        'width_median': width_median,
        'width_std': width_std,
        'width_sem': width_sem,
        'signal_mean': signal_mean,
        'signal_median': signal_median,
        'signal_std': signal_std,
        'signal_sem': signal_sem,
        'signal_col': signal_col,
        'n_cols': n_cols,
    }


# Plot combined peak-count bar chart
def plot_combined_peak_counts(df_stats):
    if df_stats.empty:
        print("No data to plot.")
        return

    df = df_stats.copy()
    df['species_order'] = df['species'].apply(lambda x: 0 if x == 'human' else 1)
    df = df.sort_values(['species_order', 'order']).reset_index(drop=True)

    x_pos   = np.arange(len(df))
    values  = df['n_peaks'].values
    labels  = df['display_name'].values

    # Assign colors by sample name
    colors = [SAMPLE_COLORS.get(name, '#BDBDBD') for name in labels]

    fig, ax = plt.subplots(figsize=(11, 6))

    # Narrow bars (0.55)
    bars = ax.bar(x_pos, values, color=colors,
                  edgecolor='black', linewidth=0.7,
                  width=0.55, zorder=3)

    # x-axis
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=10)

    # y-axis
    ax.set_ylabel('Number of detected iM peaks', fontsize=12)

    ax.set_title('Detected iM peaks across human and rice datasets',
                 fontsize=13, fontweight='bold', pad=12)

    # Light grid
    ax.grid(axis='y', linestyle='--', alpha=0.35, lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    ax.set_ylim(bottom=0)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

    # Annotate values on top of the bars
    y_max = max(values)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.012 * y_max,
                f'{int(val):,}',
                ha='center', va='bottom', fontsize=9,
                color='#333333')

    # Leave 15% headroom above the y-axis limit for value labels
    ax.set_ylim(top=y_max * 1.15)

    # Change 3: legend order — dark on top (Rice Input before Rice IgG)
    legend_elements = [
        Patch(facecolor=COLOR_HUMAN_CUTTAG,
              edgecolor='black', linewidth=0.7,
              label='Human CUT&Tag'),
        Patch(facecolor=COLOR_HUMAN_DIPSEQ,
              edgecolor='black', linewidth=0.7,
              label='Human DIP-seq'),
        Patch(facecolor=COLOR_RICE_INPUT,          # dark green
              edgecolor='black', linewidth=0.7,
              label='Rice iM-IP vs Input'),
        Patch(facecolor=COLOR_RICE_IGG,            # light green
              edgecolor='black', linewidth=0.7,
              label='Rice iM-IP vs IgG'),
    ]

    # Change 4: 4-row vertical legend with wider spacing to avoid overlap when font size changes
    ax.legend(handles=legend_elements,
              loc='upper right',          # upper right of the plot area (the rightmost bar is short, so there is room)
              # to move outside the right side of the axes, replace the line above with:
              # bbox_to_anchor=(1.02, 1.0), loc='upper left',
              ncol=1,                     # 4 vertical rows
              frameon=True,
              facecolor='white',
              edgecolor='#CCCCCC',
              framealpha=1.0,
              fontsize=9.5,
              labelspacing=0.9,           # row spacing ↑ prevents overlap
              handlelength=1.8,           # swatch length
              handleheight=1.2,           # swatch height ↑
              handletextpad=0.7,          # swatch-to-text spacing
              borderpad=0.7,              # legend padding ↑
              borderaxespad=0.6)          # legend-to-axes spacing

    plt.tight_layout()
    fn_png = os.path.join(OUTPUT_DIR, 'peak_counts_iM_combined.png')
    fn_pdf = os.path.join(OUTPUT_DIR, 'peak_counts_iM_combined.pdf')
    plt.savefig(fn_png, dpi=300, bbox_inches='tight')
    plt.savefig(fn_pdf, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Combined peak counts plot saved:")
    print(f"    {fn_png}")
    print(f"    {fn_pdf}")


# Main program
def main():
    print("===== 011_basic_stats.py =====")
    print("Scanning experiment data files...")
    files = scan_files()
    print(f"  Found {len(files)} files.")
    for f in files:
        print(f"    {f['name']}  ({f['species']}, {f['cell']})")

    print("\nProcessing files...")
    stats_list = []
    for f in files:
        print(f"  Processing {f['name']}...")
        try:
            stats = read_bed_stats(f)
            stats_list.append(stats)
            print(f"    -> {stats['n_peaks']:,} peaks")
        except Exception as e:
            print(f"    Error: {e}")

    df_stats = pd.DataFrame(stats_list)
    df_stats.to_csv(os.path.join(OUTPUT_DIR, 'basic_stats.csv'), index=False)
    print(f"\nStatistics saved to {OUTPUT_DIR}/basic_stats.csv")

    print("\nGenerating combined peak count plot...")
    plot_combined_peak_counts(df_stats)

    print(f"\nAll done! Results in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
