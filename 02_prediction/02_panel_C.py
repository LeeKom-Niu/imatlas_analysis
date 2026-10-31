#!/usr/bin/env python3
from _common import *
from matplotlib.lines import Line2D

set_style()
df = load_data()

def draw_C(ax, df, order, colors, keycol, title=''):
    x_all = df['GC_percent'].values
    y_all = df['density_Mb'].values
    rho_all, p_all = spearman_stats(x_all, y_all)

    handles = []
    labels = []

    for cat in order:
        sub = df[df[keycol] == cat]
        if len(sub) < 5: continue
        rho, p = spearman_stats(sub['GC_percent'], sub['density_Mb'])
        
        ax.scatter(sub['GC_percent'], sub['density_Mb'],
                   s=marker_size(len(sub)), alpha=0.55, color=colors[cat],
                   edgecolor='none', zorder=3)
        
        # Manually build the legend with a uniform point size of 6
        label_str = f'{LABEL_9CAT.get(cat, cat)} (r={rho:.2f}, {format_p_value(p)})'
        handle = Line2D([], [], marker='o', color='w', markerfacecolor=colors[cat], 
                        markersize=6, markeredgecolor='none', label=label_str)
        handles.append(handle)
        labels.append(label_str)

    # Fitted line
    z = np.polyfit(x_all, y_all, 1)
    xs = np.linspace(x_all.min() - 5, x_all.max() + 5, 200)
    ax.plot(xs, np.polyval(z, xs), 'k--', lw=1.4, alpha=0.85, zorder=4)
    
    all_label = f'All (r={rho_all:.2f}, {format_p_value(p_all)})'
    all_handle = Line2D([], [], color='black', linestyle='--', lw=1.4, label=all_label)
    handles.append(all_handle)
    labels.append(all_label)

    ax.set_xlabel('GC content (%)')
    ax.set_ylabel('iM density (iMs / Mb)')
    
    ax.set_ylim(-300, 2600)
    ax.set_yticks(np.arange(0, 2501, 500))
    ax.set_xlim(x_all.min() - 2, x_all.max() + 2)
    
    ax.grid(alpha=0.2, ls='--', lw=0.5); ax.set_axisbelow(True)
    ax.legend(handles=handles, labels=labels, loc='upper left', fontsize=6.5, handletextpad=0.3, frameon=False)
    
    if title: ax.set_title(title, fontsize=10, fontweight='bold', loc='left', pad=8)

# C main figure (4 domains)
fig, ax = plt.subplots(figsize=(5.5, 4.5))
draw_C(ax, df, ORDER_4, COLORS_DOMAIN, 'domain4', 'C  GC content vs iM density')
save(fig, f'{OUT_MAIN}/panel_C_4domain.png')

# C supplementary figure (9 categories)
fig, ax = plt.subplots(figsize=(6.5, 5))
draw_C(ax, df, ORDER_9, COLORS_9CAT, 'Classification_9cat', 'C  GC content vs iM density (9 cat)')
save(fig, f'{OUT_SUPP}/panel_C_9cat.png')

# C supplementary figure (6 eukaryotic subcategories)
fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
for ax, cat in zip(axes.flat, EUK_6):
    draw_C(ax, df[df['Classification_9cat'] == cat], [cat], COLORS_9CAT, 'Classification_9cat', 
           title=LABEL_9CAT.get(cat, cat)) # titles consistently use LABEL_9CAT
plt.tight_layout()
save(fig, f'{OUT_SUPP}/panel_C_euk6.png')

print('[OK] C saved')
