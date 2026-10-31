import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

# Set style (suppress palette warnings)
sns.set_theme(style="ticks", rc={
    'font.family': 'DejaVu Sans', 
    'pdf.fonttype': 42, 
    'ps.fonttype': 42
})

base = '/datapool/home/2023200496/niulk/my_project/09_type'
colors = {'Archaea & Bacteria': '#5C8F50', 'Eukaryotes': '#B65050'}

# 1. Read and process data
df_prop = pd.read_csv(f'{base}/tree_relative_proportion_filtered.tsv', sep='\t')
df_den = pd.read_csv(f'{base}/tree_total_density_filtered.tsv', sep='\t')

# Process proportion data
df_prop['Type4_5'] = (df_prop['Type 4'] + df_prop['Type 5']) * 100
df_prop['Group'] = df_prop['Domain'].apply(lambda x: 'Eukaryotes' if x == 'Eukaryote' else 'Archaea & Bacteria')

# Process density data
df_den['Group'] = df_den['Domain'].apply(lambda x: 'Eukaryotes' if x == 'Eukaryote' else 'Archaea & Bacteria')

# Figure 1: Type 4+5 proportion box plot
print("=== 1. Type 4+5 Proportion Statistics ===")
prok_prop = df_prop[df_prop['Group'] == 'Archaea & Bacteria']['Type4_5']
euk_prop  = df_prop[df_prop['Group'] == 'Eukaryotes']['Type4_5']

# Statistics and test
stat1, p1 = mannwhitneyu(prok_prop, euk_prop, alternative='two-sided')
print(f"Archaea & Bacteria (n={len(prok_prop)}): Median={prok_prop.median():.2f}%, IQR=({prok_prop.quantile(0.25):.2f}-{prok_prop.quantile(0.75):.2f})%")
print(f"Eukaryotes (n={len(euk_prop)}): Median={euk_prop.median():.2f}%, IQR=({euk_prop.quantile(0.25):.2f}-{euk_prop.quantile(0.75):.2f})%")
print(f"Mann-Whitney U test: P = {p1:.2e}\n")

fig, ax = plt.subplots(figsize=(5, 5))
sns.boxplot(data=df_prop, x='Group', y='Type4_5', hue='Group', palette=colors, 
            width=0.5, showfliers=False, legend=False, ax=ax)
sns.stripplot(data=df_prop, x='Group', y='Type4_5', 
              color='black', alpha=0.6, size=4, jitter=0.15, ax=ax)

# Significance annotation
y_max = df_prop['Type4_5'].max()
y_offset = y_max * 0.05
ax.plot([0, 0, 1, 1], [y_max + y_offset, y_max + y_offset*2, y_max + y_offset*2, y_max + y_offset], lw=1.5, c='black')
ax.text(0.5, y_max + y_offset*2.2, f'***P = {p1:.2e}', ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Proportion of Type 4+5 iMs (%)', fontsize=12, fontweight='bold')
ax.set_xlabel('')
ax.set_title('Structural complexity of iMs', pad=15, fontsize=14, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(f'{base}/Supplementary_Fig_Type4_5_Boxplot.pdf', bbox_inches='tight')
plt.savefig(f'{base}/Supplementary_Fig_Type4_5_Boxplot.png', dpi=300, bbox_inches='tight')
plt.close()


# Figure 2: Total iM density box plot
print("=== 2. Total iM Density Statistics ===")
prok_den = df_den[df_den['Group'] == 'Archaea & Bacteria']['Total_Density']
euk_den  = df_den[df_den['Group'] == 'Eukaryotes']['Total_Density']

stat2, p2 = mannwhitneyu(prok_den, euk_den, alternative='two-sided')
print(f"Archaea & Bacteria (n={len(prok_den)}): Median={prok_den.median():.2f}, IQR=({prok_den.quantile(0.25):.2f}-{prok_den.quantile(0.75):.2f})")
print(f"Eukaryotes (n={len(euk_den)}): Median={euk_den.median():.2f}, IQR=({euk_den.quantile(0.25):.2f}-{euk_den.quantile(0.75):.2f})")
print(f"Mann-Whitney U test: P = {p2:.2e}\n")

fig, ax = plt.subplots(figsize=(5, 5))
sns.boxplot(data=df_den, x='Group', y='Total_Density', hue='Group', palette=colors, 
            width=0.5, showfliers=False, legend=False, ax=ax)
sns.stripplot(data=df_den, x='Group', y='Total_Density', 
              color='black', alpha=0.6, size=4, jitter=0.15, ax=ax)

# Significance annotation
y_max_den = df_den['Total_Density'].max()
y_offset_den = y_max_den * 0.05
ax.plot([0, 0, 1, 1], [y_max_den + y_offset_den, y_max_den + y_offset_den*2, y_max_den + y_offset_den*2, y_max_den + y_offset_den], lw=1.5, c='black')
ax.text(0.5, y_max_den + y_offset_den*2.2, f'ns (P = {p2:.3f})' if p2 >= 0.05 else f'*P = {p2:.2e}', ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Total iM density (iMs/Mb)', fontsize=12, fontweight='bold')
ax.set_xlabel('')
ax.set_title('Total iM density', pad=15, fontsize=14, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(f'{base}/Supplementary_Fig_Total_Density_Boxplot.pdf', bbox_inches='tight')
plt.savefig(f'{base}/Supplementary_Fig_Total_Density_Boxplot.png', dpi=300, bbox_inches='tight')
plt.close()

print("✅ Both supplementary figures saved!")
