#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os

BASE_DIR = '/datapool/home/2023200496/niulk/my_project/09_type'
CLEAN_FILE = '/datapool/home/2023200496/niulk/my_project/20_iM_DB_analysis/02_correlation_plots/data/iM_summary_clean.tsv'

# 1. Read iM count matrices for the three domains of life
dfs = []
for folder in ['01_archaea', '02_bacteria', '03_eukaryote']:
    file_path = os.path.join(BASE_DIR, folder, 'imotif_counts_matrix.tsv')
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t')
        df['Domain'] = folder.split('_')[1].capitalize()  # Archaea, Bacteria, Eukaryote
        dfs.append(df)
        print(f"Loaded {len(df)} species from {folder}")
    else:
        print(f"[WARN] File not found: {file_path}")

counts = pd.concat(dfs, ignore_index=True)
types = ['Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5']

# 2. Read genome length information
clean = pd.read_csv(CLEAN_FILE, sep='\t')[['Genome', 'genome_length_bp']]

# 3. Merge and compute
merged = counts.merge(clean, on='Genome', how='left')
missing_len = merged['genome_length_bp'].isna().sum()
if missing_len > 0:
    print(f"[WARNING] {missing_len} 个物种在 clean 表中找不到基因组长度，将导致绝对密度为 NaN！")
    print(merged[merged['genome_length_bp'].isna()]['Genome'].tolist()[:10], "...")

# Compute totals
merged['Total_Count'] = merged[types].sum(axis=1)

# Generate dataset 1: 100% relative proportions (for stacked bar chart)
rel_data = merged[['Genome', 'Domain']].copy()
for t in types:
    rel_data[t] = merged[t] / merged['Total_Count']
# Handle division by zero (species without iM)
rel_data[types] = rel_data[types].fillna(0)
rel_out = os.path.join(BASE_DIR, 'tree_relative_proportion.tsv')
rel_data.to_csv(rel_out, sep='\t', index=False)
print(f"[OK] 相对比例数据已保存: {rel_out}")

# Generate dataset 2: absolute density iMs/Mb (alternative for stacked bar chart)
abs_data = merged[['Genome', 'Domain']].copy()
for t in types:
    # Divide by 2 (deduplicate across both strands), then by genome length (Mb)
    abs_data[t] = (merged[t] / 2) / (merged['genome_length_bp'] / 1e6)
abs_out = os.path.join(BASE_DIR, 'tree_absolute_density.tsv')
abs_data.to_csv(abs_out, sep='\t', index=False)
print(f"[OK] 绝对密度数据已保存: {abs_out}")

# Generate dataset 3: total density log10 (for heatmap/bubble plot on the outer ring of the tree)
total_data = abs_data[['Genome', 'Domain']].copy()
total_data['Total_Density'] = abs_data[types].sum(axis=1)
total_data['Total_log10'] = np.log10(total_data['Total_Density'] + 1)
total_out = os.path.join(BASE_DIR, 'tree_total_density.tsv')
total_data[['Genome', 'Domain', 'Total_Density', 'Total_log10']].to_csv(total_out, sep='\t', index=False)
print(f"[OK] 总密度/热图数据已保存: {total_out}")
