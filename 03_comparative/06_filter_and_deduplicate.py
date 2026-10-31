#!/usr/bin/env python3
import pandas as pd
import os

BASE_DIR = '/datapool/home/2023200496/niulk/my_project/09_type'

# 1. Read total density file
print("正在读取总密度文件进行去重和名称清洗...")
df_total = pd.read_csv(f'{BASE_DIR}/tree_total_density.tsv', sep='\t')

# Filter out entries with zero density
df_total = df_total[df_total['Total_Density'] > 0].copy()

# Extract clean Genus_species name (first two words)
df_total['CleanGenome'] = df_total['Genome'].apply(lambda x: '_'.join(x.split('_')[:2]))

# 2. Sort: by density from high to low (so deduplication keeps the highest-density version)
df_total = df_total.sort_values(by='Total_Density', ascending=False)

# 3. Special-case human: if 'Homo_sapiens' exists, force-prioritize GRCh38_p14
human_mask = df_total['CleanGenome'] == 'Homo_sapiens'
if human_mask.any():
    # Find the GRCh38 row
    human38_mask = df_total['Genome'] == 'Homo_sapiens_GRCh38_p14'
    if human38_mask.any():
        # Move GRCh38 to the front to ensure it is retained during deduplication
        human38 = df_total[human38_mask]
        others = df_total[~human38_mask]
        df_total = pd.concat([human38, others]).reset_index(drop=True)

# 4. Core deduplication: keep the first per CleanGenome (highest density or the specially approved GRCh38)
df_dedup = df_total.drop_duplicates(subset=['CleanGenome'], keep='first').copy()

# Build dict: original long name -> clean short name
name_map = dict(zip(df_dedup['Genome'], df_dedup['CleanGenome']))
keep_genomes = set(df_dedup['Genome'])

print(f"去重完成！原始物种数: {len(df_total)}，去重后保留: {len(df_dedup)} 个")

# 5. Apply cleaning and filtering to all three files
files_to_filter = [
    'tree_relative_proportion.tsv',
    'tree_absolute_density.tsv',
    'tree_total_density.tsv'
]

for fname in files_to_filter:
    filepath = f'{BASE_DIR}/{fname}'
    df = pd.read_csv(filepath, sep='\t')
    
    # Filter: keep only retained genomes
    df_filtered = df[df['Genome'].isin(keep_genomes)].copy()
    
    # Rename: replace long names with clean Genus_species names
    df_filtered['Genome'] = df_filtered['Genome'].map(name_map)
    
    # Re-sort: first by Domain, then by clean species name
    df_filtered = df_filtered.sort_values(by=['Domain', 'Genome']).reset_index(drop=True)
    
    # Save
    out_path = f'{BASE_DIR}/{fname.replace(".tsv", "_filtered.tsv")}'
    df_filtered.to_csv(out_path, sep='\t', index=False)
    print(f"已保存: {out_path} (行数: {len(df_filtered)})")

# 6. Also output the list of retained clean names to help check the tree file
df_dedup[['CleanGenome', 'Domain']].to_csv(f'{BASE_DIR}/kept_species_names.txt', sep='\t', index=False, header=False)
print("\n全部完成！请使用带 _filtered 后缀的三个文件去网页工具画图。")
print("同时输出了 kept_species_names.txt，你可以用它核对树文件的叶子名。")
