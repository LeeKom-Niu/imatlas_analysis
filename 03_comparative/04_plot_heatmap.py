#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

sys.path.append('/datapool/home/2023200496/niulk/my_project/20_iM_DB_analysis/02_correlation_plots')
from _common import set_style, OUT_SUPP

set_style(fontsize=10)

IN_FILE = '../results/all_species_region_density.tsv'
OUT_DIR = '../results'
df = pd.read_csv(IN_FILE, sep='\t')
print(f'Loaded {len(df)} species')

# 1. Compute Fold Enrichment
regions = ['Prom', 'Gene', 'Down', 'Distal']
for r in regions:
    # fold = region density / whole-genome density
    df[f'{r}_fold'] = df[f'{r}_density'] / df['Genome_density']
    
# 2. Join 9cat classification and filter to eukaryotes
CLEAN = '/datapool/home/2023200496/niulk/my_project/20_iM_DB_analysis/02_correlation_plots/data/iM_summary_clean.tsv'
clean_df = pd.read_csv(CLEAN, sep='\t')[['Genome', 'Classification_9cat']]
df = df.merge(clean_df, left_on='Species', right_on='Genome', how='left')

# Keep only the 6 eukaryotic categories
EUK_6 = ['Protozoa', 'Fungi', 'Plant', 'Invertebrate', 'Vertebrate Other', 'Mammalian']
df = df[df['Classification_9cat'].isin(EUK_6)]

print(f'After mapping and filtering: {len(df)} species')

# New: save raw Fold Enrichment data for each species
species_fold_out = f'{OUT_DIR}/all_species_fold_enrichment.tsv'
df[['Species', 'Classification_9cat', 'Genome_density', 
    'Prom_fold', 'Gene_fold', 'Down_fold', 'Distal_fold']].to_csv(species_fold_out, sep='\t', index=False)
print(f'[OK] 每个物种的 Fold 数据已保存至: {species_fold_out}')

# 3. Compute summary statistics per category and region (Median, Mean, Q1, Q3)
summary = []
for cat in EUK_6:
    cat_df = df[df['Classification_9cat'] == cat]
    if len(cat_df) == 0: continue
    
    row = {'Category': cat, 'n': len(cat_df)}
    for r in regions:
        row[f'{r}_median'] = cat_df[f'{r}_fold'].median()
        row[f'{r}_mean'] = cat_df[f'{r}_fold'].mean()
        row[f'{r}_q1'] = cat_df[f'{r}_fold'].quantile(0.25)
        row[f'{r}_q3'] = cat_df[f'{r}_fold'].quantile(0.75)
    summary.append(row)

sum_df = pd.DataFrame(summary).set_index('Category').reindex(EUK_6)

# New: save the summary table behind the heatmap
summary_out = f'{OUT_DIR}/heatmap_summary_stats.tsv'
sum_df.to_csv(summary_out, sep='\t')
print(f'[OK] 热图汇总统计（Median/Mean/IQR）已保存至: {summary_out}')

# 4. Draw heatmap (using median data)
plot_data = sum_df[[f'{r}_median' for r in regions]].copy()
plot_data.columns = ['Promoter', 'Gene body', 'Downstream', 'Distal']

fig, ax = plt.subplots(figsize=(6, 4.5))
log_data = np.log2(plot_data)

sns.heatmap(log_data, 
            cmap='RdBu_r', center=0, 
            annot=plot_data, fmt='.2f', # annotate raw fold values
            cbar_kws={'label': 'log2(Fold enrichment)'},
            linewidths=0.5, linecolor='white',
            ax=ax)

ax.set_title('iM enrichment across genomic features', loc='left', pad=12)
ax.set_ylabel('')
ax.set_xlabel('')

plt.tight_layout()
out_path = f'{OUT_SUPP}/panel_heatmap_enrichment.png'
fig.savefig(out_path, dpi=300, bbox_inches='tight')
fig.savefig(out_path.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
print(f'[OK] 热图已保存至: {out_path}')
