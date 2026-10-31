#!/usr/bin/env python3
"""
viral_i-motif_density_by_host.py  (v2 - Nature style + narrow bars + iM)

Viral iM density by host type (eukaryote vs prokaryote)

Output:
  viral_iM_density_by_host_euk_prok.png / .pdf

Font: pdf.fonttype = 42 (editable in Illustrator)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from scipy.stats import mannwhitneyu


# Nature global style + editable PDF
mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.labelsize': 11.5,
    'axes.titlesize': 12,
    'axes.linewidth': 0.9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.labelsize': 10.5,
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


# Color scheme (low-saturation blue vs green)
HOST_COLORS = {
    'eukaryote':   '#3A6B85',   # dark blue
    'prokaryote':  '#7BA05B',   # olive green
}

HOST_LABELS = {
    'eukaryote':   'Eukaryotic host',
    'prokaryote':  'Prokaryotic host',
}


# Read data
file_path = 'viral_with_host.tsv'
data = pd.read_csv(file_path, sep='\t')

target_hosts = ['eukaryote', 'prokaryote']
filtered_data = data[data['Host'].isin(target_hosts)].copy()

print("Data preview:")
print(filtered_data.head())
print("\nHost distribution:")
print(filtered_data['Host'].value_counts())


# Statistics
host_stats = (
    filtered_data
    .groupby('Host')['Genomic_density']
    .agg(mean='mean', std='std', count='count')
    .reset_index()
)
host_stats['sem'] = host_stats['std'] / np.sqrt(host_stats['count'])
host_stats = host_stats.set_index('Host').reindex(target_hosts).reset_index()

print("\nHost statistics:")
print(host_stats.to_string(index=False))


# Significance test
groups = {}
for host in target_hosts:
    groups[host] = filtered_data[
        filtered_data['Host'] == host
    ]['Genomic_density'].values

stat_u, p_value = mannwhitneyu(groups['eukaryote'],
                                groups['prokaryote'],
                                alternative='two-sided')

if p_value < 0.001:
    sig = '***'
elif p_value < 0.01:
    sig = '**'
elif p_value < 0.05:
    sig = '*'
else:
    sig = 'ns'

print(f"\nMann-Whitney U: U = {stat_u:.1f}, P = {p_value:.3e}  ({sig})")


# Plotting
fig, ax = plt.subplots(figsize=(5.0, 5.2))

x_pos = np.arange(len(target_hosts))
means = host_stats['mean'].values
sems  = host_stats['sem'].values
colors = [HOST_COLORS[h] for h in target_hosts]

# Narrow bars (width=0.42)
bars = ax.bar(x_pos, means,
               yerr=sems,
               width=0.42,
               color=colors,
               edgecolor='black', linewidth=0.8,
               capsize=5,
               error_kw=dict(elinewidth=1.2, capthick=1.2,
                             ecolor='#333333'),
               zorder=3)

# x-axis labels (with n)
xtick_labels = []
for h in target_hosts:
    n = int(host_stats.loc[host_stats['Host'] == h, 'count'].values[0])
    xtick_labels.append(f'{HOST_LABELS[h]}\n(n={n:,})')

ax.set_xticks(x_pos)
ax.set_xticklabels(xtick_labels, fontsize=10.5)

# y-axis
ax.set_ylabel('iM density (per Mb)', fontsize=11.5)

# y-axis upper limit: reserve 40% headroom for the significance marker
y_max_data = (means + sems).max()
ax.set_ylim(0, y_max_data * 1.40)

# Grid
ax.grid(axis='y', ls='--', lw=0.6, alpha=0.25, zorder=0)
ax.set_axisbelow(True)

# Significance line + asterisk
y_line = y_max_data * 1.15
y_text = y_max_data * 1.20
ax.plot([x_pos[0], x_pos[0], x_pos[1], x_pos[1]],
        [y_line, y_line + y_max_data * 0.02,
         y_line + y_max_data * 0.02, y_line],
        color='black', lw=1.0, zorder=5)
ax.text((x_pos[0] + x_pos[1]) / 2, y_text, sig,
        ha='center', va='bottom',
        fontsize=13, fontweight='bold', zorder=5)

# Title
ax.set_title('Viral iM density by host type',
             fontsize=12, fontweight='bold', loc='left', pad=10)

plt.tight_layout()
fig.savefig('viral_iM_density_by_host_euk_prok.png',
            dpi=300, bbox_inches='tight')
fig.savefig('viral_iM_density_by_host_euk_prok.pdf',
            dpi=300, bbox_inches='tight')
plt.close(fig)

print("\nsaved: viral_iM_density_by_host_euk_prok.png / .pdf")


# Descriptive statistics
print("\n" + "=" * 55)
print("Descriptive statistics")
print("=" * 55)
for h in target_hosts:
    d = filtered_data[filtered_data['Host'] == h]['Genomic_density']
    print(f"\n{h}:")
    print(f"  n      : {len(d):,}")
    print(f"  mean   : {d.mean():.2f}")
    print(f"  SD     : {d.std():.2f}")
    print(f"  SEM    : {d.sem():.2f}")
    print(f"  median : {d.median():.2f}")
    print(f"  range  : {d.min():.2f} – {d.max():.2f}")
