import pandas as pd
from scipy.stats import mannwhitneyu

# 1. Read data
base = '/datapool/home/2023200496/niulk/my_project/09_type'
df_prop = pd.read_csv(f'{base}/tree_relative_proportion_filtered.tsv', sep='\t')
df_den = pd.read_csv(f'{base}/tree_total_density_filtered.tsv', sep='\t')

types = ['Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5']

# Part 1: Overall mean proportion of iM types (98 species mean)
print("=== 1. Overall Proportions (98 species mean) ===")
overall_prop = (df_prop[types].mean() * 100).round(2)
print(overall_prop)

# Part 2: Type 4+5 proportion (Prokaryotes vs Eukaryotes)
df_prop['Type4_5'] = df_prop['Type 4'] + df_prop['Type 5']

# Grouping: prokaryotes (Archaea + Bacteria) vs eukaryotes (Eukaryote)
prok = df_prop[df_prop['Domain'].isin(['Archaea', 'Bacteria'])]['Type4_5']
euk  = df_prop[df_prop['Domain'] == 'Eukaryote']['Type4_5']

print("\n=== 2. Type 4+5 Proportion ===")
print(f"Prokaryotes (Archaea+Bacteria, n={len(prok)}):")
print(f"  Mean: {prok.mean()*100:.2f}%")
print(f"  Median: {prok.median()*100:.2f}%")
print(f"  IQR: {prok.quantile(0.25)*100:.2f}% - {prok.quantile(0.75)*100:.2f}%")

print(f"Eukaryotes (n={len(euk)}):")
print(f"  Mean: {euk.mean()*100:.2f}%")
print(f"  Median: {euk.median()*100:.2f}%")
print(f"  IQR: {euk.quantile(0.25)*100:.2f}% - {euk.quantile(0.75)*100:.2f}%")

# Mann-Whitney U test (two-sided, the most rigorous)
stat, p_value = mannwhitneyu(prok, euk, alternative='two-sided')
print(f"\nMann-Whitney U test (two-sided): U={stat:.1f}, P = {p_value:.2e}")

# Part 3: Total density extremes
print("\n=== 3. Total density extremes ===")
min_row = df_den.loc[df_den['Total_Density'].idxmin()]
max_row = df_den.loc[df_den['Total_Density'].idxmax()]
print(f"Min: {min_row['Genome']} ({min_row['Domain']}) = {min_row['Total_Density']:.2f}")
print(f"Max: {max_row['Genome']} ({max_row['Domain']}) = {max_row['Total_Density']:.2f}")

# Part 4: Median total density, prokaryotes vs eukaryotes (inner-ring description of supplementary figure 4D)
print("\n=== 4. Total density median (Prok vs Euk) ===")
den_prok = df_den[df_den['Domain'].isin(['Archaea', 'Bacteria'])]['Total_Density']
den_euk  = df_den[df_den['Domain'] == 'Eukaryote']['Total_Density']
print(f"Prokaryotes median density: {den_prok.median():.2f} iMs/Mb")
print(f"Eukaryotes median density: {den_euk.median():.2f} iMs/Mb")
