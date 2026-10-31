#!/usr/bin/env python3
"""
03_genomic_region_stacked_bar.py  (v3 - main and supplementary figures output separately, sample names renamed)

Four-region classification of iM peaks per dataset (mutually exclusive, summing to 100%):
  1. Promoter          : 1 kb upstream of TSS
  2. Gene body         : TSS to TES
  3. Downstream        : 1 kb downstream of TES
  4. Distal intergenic : outside the above

Priority: Promoter > Gene body > Downstream > Intergenic

Output:
  Main figure: datasets specified by --main
  Supplementary figure: datasets specified by --supp

Usage:
  cd ~/niulk/my_project/19_modified_analysis/10_iMab/05_Experiment_analysis

  python3 03_genomic_region_stacked_bar.py \\
      --genome_hg38 ./genome \\
      --genome_hg19 ./genome_hg19 \\
      --genome_rice ./genome \\
      --main HEK293T_CUTTag hg38 ./experiment_data/human_cutntag/HEK293T_iM_highConf_intersect.bed \\
      --main WDLPS_CUTTag   hg38 ./experiment_data/human_cutntag/WDLPS_iM_highConf_intersect.bed \\
      --main Rice_IgG       rice ./experiment_data/rice/IP_vs_IgG_intersect.narrowPeak \\
      --main Rice_Input     rice ./experiment_data/rice/IP_vs_Input_intersect.narrowPeak \\
      --supp HEK293T_DIPseq hg19 ./experiment_data/human/GSE248445_HEK293T_1_2_intersect_peaks.bed.gz \\
      --supp MCF7_DIPseq    hg19 ./experiment_data/human/GSE248445_MCF7_1_2_intersect_peaks.bed.gz \\
      --supp U2OS_DIPseq    hg19 ./experiment_data/human/GSE248445_U2OS_1_2_intersect_peaks.bed.gz \\
      --supp 3cell_DIPseq   hg19 ./experiment_data/human/GSE248445_iMab_MCF7_HEK293T_U2OS_intersect_peaks.bed.gz \\
      --outdir ./03_regional_stackbar
"""

import os
import gzip
import bisect
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl


# Sample display-name mapping (command-line name → display name on plot)
DISPLAY_NAMES = {
    'HEK293T_CUTTag':  'HEK293T CUT&Tag',
    'WDLPS_CUTTag':    'WDLPS CUT&Tag',
    'HEK293T_DIPseq':  'HEK293T DIP-seq',
    'MCF7_DIPseq':     'MCF7 DIP-seq',
    'U2OS_DIPseq':     'U2OS DIP-seq',
    '3cell_DIPseq':    '3-cell DIP-seq',
    'Rice_IgG':        'Rice_IgG (pH7.0)',
    'Rice_Input':      'Rice_Input (seedling)',
}


def display_name(name):
    return DISPLAY_NAMES.get(name, name)


# IO
def open_maybe_gz(path):
    path = str(path)
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_bed3(path):
    rows = []
    with open_maybe_gz(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            try:
                rows.append((parts[0], int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    if not rows:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])


def read_bed6(path):
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


# Interval lookup
class IntervalLookup:
    __slots__ = ('starts', 'prefix_max_end')

    def __init__(self, intervals):
        self.starts = [s for s, _ in intervals]
        pme, mx = [], -1
        for _, e in intervals:
            if e > mx:
                mx = e
            pme.append(mx)
        self.prefix_max_end = pme

    def overlaps(self, qs, qe):
        if not self.starts:
            return False
        i = bisect.bisect_left(self.starts, qe)
        if i == 0:
            return False
        return self.prefix_max_end[i - 1] > qs


def build_lookups(df):
    out = {}
    for chrom, grp in df.groupby('chrom', sort=False):
        ivs = sorted(zip(grp['start'].tolist(), grp['end'].tolist()))
        out[chrom] = IntervalLookup(ivs)
    return out


def merge_intervals(df):
    if df.empty:
        return df
    rows = []
    for chrom, grp in df.groupby('chrom', sort=False):
        ivs = sorted(zip(grp['start'].tolist(), grp['end'].tolist()))
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                rows.append((chrom, cs, ce))
                cs, ce = s, e
        rows.append((chrom, cs, ce))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])


# Build 3 sets of intervals from genes
def build_regions_from_genes(genes_df, chrom_sizes,
                             promoter_len=1000,
                             downstream_len=1000):
    prom_rows = []
    gene_rows = []
    down_rows = []

    for chrom, gs, ge, strand in zip(
            genes_df['chrom'], genes_df['start'],
            genes_df['end'], genes_df['strand']):

        clen = chrom_sizes.get(chrom, 0)
        if clen <= 0:
            continue

        if strand == '+':
            tss = gs
            tes = ge
            ps = max(0, tss - promoter_len)
            pe = min(clen, tss)
            ds = max(0, tes)
            de = min(clen, tes + downstream_len)
        elif strand == '-':
            tss = ge
            tes = gs
            ps = max(0, tss)
            pe = min(clen, tss + promoter_len)
            ds = max(0, tes - downstream_len)
            de = min(clen, tes)
        else:
            continue

        prom_rows.append((chrom, ps, pe))
        gene_rows.append((chrom, gs, ge))
        down_rows.append((chrom, ds, de))

    prom_df = pd.DataFrame(prom_rows, columns=['chrom', 'start', 'end'])
    gene_df = pd.DataFrame(gene_rows, columns=['chrom', 'start', 'end'])
    down_df = pd.DataFrame(down_rows, columns=['chrom', 'start', 'end'])
    return prom_df, gene_df, down_df


# Four-region classification (mutually exclusive)
def classify_peaks(peak_df, prom_lk, gene_lk, down_lk):
    n_prom = n_gene = n_down = n_inter = 0
    for chrom, s, e in zip(peak_df['chrom'], peak_df['start'], peak_df['end']):
        lk = prom_lk.get(chrom)
        if lk is not None and lk.overlaps(s, e):
            n_prom += 1
            continue
        lk = gene_lk.get(chrom)
        if lk is not None and lk.overlaps(s, e):
            n_gene += 1
            continue
        lk = down_lk.get(chrom)
        if lk is not None and lk.overlaps(s, e):
            n_down += 1
            continue
        n_inter += 1
    return n_prom, n_gene, n_down, n_inter


# Build region cache
def build_region_cache(genome_dirs, species_map, promoter_len, downstream_len):
    cache = {}
    for coord, gdir in genome_dirs.items():
        if gdir is None:
            continue
        species = species_map[coord]
        gd = Path(gdir)
        chrom_sizes = read_chrom_sizes(gd / f'{species}.chrom.sizes')
        genes_bed = read_bed6(gd / f'{species}_genes.bed')

        if not chrom_sizes or genes_bed.empty:
            print(f'  [warn] {coord} missing chrom.sizes or genes.bed, skipping')
            continue

        prom_df, gene_df, down_df = build_regions_from_genes(
            genes_bed, chrom_sizes,
            promoter_len=promoter_len,
            downstream_len=downstream_len)

        prom_m = merge_intervals(prom_df)
        gene_m = merge_intervals(gene_df)
        down_m = merge_intervals(down_df)

        cache[coord] = {
            'prom_lk': build_lookups(prom_m),
            'gene_lk': build_lookups(gene_m),
            'down_lk': build_lookups(down_m),
        }
        print(f'[{coord}] regions built: '
              f'prom={len(prom_m)}  gene={len(gene_m)}  down={len(down_m)}')
    return cache


# Classify one dataset
def process_one(name, coord, peak_path, region_cache):
    if coord not in region_cache:
        print(f'  [skip] {name}: no genome directory provided for coordinate {coord}')
        return None

    peaks = read_bed3(peak_path)
    if peaks.empty:
        print(f'  [skip] {name}: empty file')
        return None

    r = region_cache[coord]
    n_prom, n_gene, n_down, n_inter = classify_peaks(
        peaks, r['prom_lk'], r['gene_lk'], r['down_lk'])

    total = len(peaks)
    row = {
        'sample':         display_name(name),   # ← use display name
        'sample_raw':     name,                  # keep the raw name
        'coord':          coord,
        'n_peaks':        total,
        'n_promoter':     n_prom,
        'n_gene':         n_gene,
        'n_downstream':   n_down,
        'n_intergenic':   n_inter,
        'pct_promoter':   n_prom / total * 100,
        'pct_gene':       n_gene / total * 100,
        'pct_downstream': n_down / total * 100,
        'pct_intergenic': n_inter / total * 100,
    }
    print(f'  {name:35s}  n={total:>7d}  '
          f'prom={n_prom/total*100:5.1f}%  '
          f'gene={n_gene/total*100:5.1f}%  '
          f'down={n_down/total*100:5.1f}%  '
          f'inter={n_inter/total*100:5.1f}%')
    return row


# Draw stacked bar chart
def plot_stacked(df, outdir, filename, fig_width=9,
                 show_legend=True, title=None):
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.linewidth': 0.9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'pdf.compression': 6,
    })

    df_plot = df.iloc[::-1].reset_index(drop=True)
    samples = df_plot['sample'].tolist()

    region_names = ['Promoter', 'Gene body',
                    'Downstream', 'Distal intergenic']
    region_cols = ['pct_promoter', 'pct_gene',
                   'pct_downstream', 'pct_intergenic']
    colors = ['#C44E52', '#4C72B0', '#55A868', '#BDBDBD']

    fig, ax = plt.subplots(figsize=(fig_width,
                                    max(3.0, 0.55 * len(samples) + 1.8)))

    y_pos = np.arange(len(samples))
    left = np.zeros(len(samples))

    for col, color, rname in zip(region_cols, colors, region_names):
        vals = df_plot[col].values
        ax.barh(y_pos, vals, left=left, color=color,
                edgecolor='white', linewidth=0.7,
                label=rname, height=0.7)
        for i, (v, l) in enumerate(zip(vals, left)):
            if v >= 5:
                ax.text(l + v / 2, i, f'{v:.1f}',
                        ha='center', va='center',
                        fontsize=8.5, color='white', fontweight='bold')
        left = left + vals

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples, fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_xlabel('Percentage of peaks (%)', fontsize=11)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.grid(axis='x', ls='--', lw=0.6, alpha=0.3)
    ax.set_axisbelow(True)

    if show_legend:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10),
                  ncol=4, frameon=False, fontsize=9.5)

    if title:
        ax.set_title(title, fontsize=11.5, loc='left', pad=10)

    plt.tight_layout()
    png = Path(outdir) / f'{filename}.png'
    pdf = Path(outdir) / f'{filename}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'saved: {png} / {pdf}')


# main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genome_hg38', default=None)
    ap.add_argument('--genome_hg19', default=None)
    ap.add_argument('--genome_rice', default=None)
    ap.add_argument('--main', nargs=3, action='append', default=[],
                    metavar=('NAME', 'COORD', 'PEAK_BED'),
                    help='main figure datasets')
    ap.add_argument('--supp', nargs=3, action='append', default=[],
                    metavar=('NAME', 'COORD', 'PEAK_BED'),
                    help='supplementary figure datasets')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--promoter_len', type=int, default=1000)
    ap.add_argument('--downstream_len', type=int, default=1000)
    args = ap.parse_args()

    if not args.main and not args.supp:
        raise SystemExit('At least one of --main or --supp must be specified')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Region cache
    genome_dirs = {
        'hg38': args.genome_hg38,
        'hg19': args.genome_hg19,
        'rice': args.genome_rice,
    }
    species_map = {'hg38': 'human', 'hg19': 'human', 'rice': 'rice'}

    print('=' * 60)
    print('Building region annotations...')
    region_cache = build_region_cache(
        genome_dirs, species_map,
        args.promoter_len, args.downstream_len)

    # Main figure
    if args.main:
        print('\n' + '=' * 60)
        print('Processing MAIN figure datasets...')
        rows_main = []
        for name, coord, path in args.main:
            r = process_one(name, coord, path, region_cache)
            if r is not None:
                rows_main.append(r)

        if rows_main:
            df_main = pd.DataFrame(rows_main)
            df_main.to_csv(outdir / 'region_stats_main.tsv',
                           sep='\t', index=False)
            plot_stacked(df_main, outdir, 'region_stacked_bar_main',
                         fig_width=8, show_legend=True)
            print(f'\nsaved: {outdir}/region_stats_main.tsv')

    # Supplementary figure
    if args.supp:
        print('\n' + '=' * 60)
        print('Processing SUPPLEMENTARY figure datasets...')
        rows_supp = []
        for name, coord, path in args.supp:
            r = process_one(name, coord, path, region_cache)
            if r is not None:
                rows_supp.append(r)

        if rows_supp:
            df_supp = pd.DataFrame(rows_supp)
            df_supp.to_csv(outdir / 'region_stats_supp.tsv',
                           sep='\t', index=False)
            plot_stacked(df_supp, outdir, 'region_stacked_bar_supp',
                         fig_width=8, show_legend=True)
            print(f'\nsaved: {outdir}/region_stats_supp.tsv')

    print('\nAll done!')


if __name__ == '__main__':
    main()
