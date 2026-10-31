#!/usr/bin/env python3
"""
12_density_barplot_9cat.py
Generate a 9-category density horizontal bar plot from the latest curated data (iM_summary_clean.tsv).
Includes:
- Horizontal bar plot (sorted by mean density, descending)
- Single-hue gradient color scheme
- Mean ± SEM error bars
- Value labels
"""
from _common import *
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

set_style(fontsize=10)
df = load_data()

# Statistics: all genomes (compute mean, std, SEM)
cat_stats = df.groupby('Classification_9cat')['density_Mb'].agg(
    n='count',
    mean='mean',
    std='std',
    n_zero=lambda x: (x == 0).sum(),
).reset_index()

cat_stats['sem'] = cat_stats['std'] / np.sqrt(cat_stats['n'])
cat_stats['pct_zero'] = cat_stats['n_zero'] / cat_stats['n'] * 100

# Sort by mean density in descending order
cat_stats = cat_stats.sort_values('mean', ascending=False).reset_index(drop=True)

print('=== Density statistics (ALL genomes, sorted by mean) ===')
print(cat_stats[['Classification_9cat', 'n', 'mean', 'std', 'sem', 'pct_zero']]
      .to_string(index=False))

# Color scheme: single-hue gradient (dark blue → light blue)
BLUES_CMAP = LinearSegmentedColormap.from_list(
    'blue_gradient',
    ['#1F4E79', '#2E75B6', '#5B9BD5', '#9DC3E6', '#D6E4F0']
)
n_cats = len(cat_stats)
bar_colors = [BLUES_CMAP(i / max(n_cats - 1, 1)) for i in range(n_cats)]

# Plotting
fig, ax = plt.subplots(figsize=(8.5, 5))

y_pos = np.arange(n_cats)[::-1]
means  = cat_stats['mean'].values
sems   = cat_stats['sem'].values
labels = [LABEL_9CAT.get(c, c) for c in cat_stats['Classification_9cat']]

# Bar chart + SEM
bars = ax.barh(y_pos, means,
               xerr=sems,
               color=bar_colors,
               height=0.6,
               edgecolor='#333333', linewidth=0.8,
               capsize=3,
               error_kw=dict(elinewidth=1.0, capthick=1.0, ecolor='#333333'),
               zorder=3)

# Value labels (to the right of error bars)
for bar, mean_v, sem_v in zip(bars, means, sems):
    x_text = mean_v + sem_v + means.max() * 0.025
    ax.text(x_text,
            bar.get_y() + bar.get_height() / 2,
            f'{mean_v:.2f}',
            ha='left', va='center',
            fontsize=9.5, color='#222222')

# Axis settings
ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xlabel('Mean iM density (iMs / Mb)')

x_max = (means + sems).max()
ax.set_xlim(0, x_max * 1.22)

# Grid and title
ax.grid(axis='x', ls='--', lw=0.5, alpha=0.25, zorder=0)
ax.set_axisbelow(True)
ax.set_title('Mean iM density by category\n(all genomes; mean ± SEM)',
             loc='left', pad=10)

plt.tight_layout()

# Save output
save(fig, f'{OUT_SUPP}/panel_density_9cat.png')

# Save statistics table
out_table = cat_stats.rename(columns={'Classification_9cat': 'Category'})
out_table.to_csv(f'{OUT_SUPP}/density_stats_9cat.tsv', sep='\t', index=False)

print(f'\n[OK] Saved: {OUT_SUPP}/panel_density_9cat.png')
print(f'[OK] Saved: {OUT_SUPP}/density_stats_9cat.tsv')
