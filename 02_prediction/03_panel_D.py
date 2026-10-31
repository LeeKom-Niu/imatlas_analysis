#!/usr/bin/env python3
from _common import *
from matplotlib.lines import Line2D

set_style()
df = load_data()

def draw_D(ax, df, order, colors, keycol, title=''):
    x_all = df['log10_size'].values
    y_all = df['log10_im1'].values

    handles = []
    labels = []

    for cat in order:
        sub = df[df[keycol] == cat]
        if len(sub) < 5: continue
        rho, p = spearman_stats(sub['log10_size'], sub['log10_im1'])
        
        ax.scatter(sub['log10_size'], sub['log10_im1'],
                   s=marker_size(len(sub)), alpha=0.55, color=colors[cat],
                   edgecolor='none', zorder=3)
        
        # Manually build the legend with a uniform point size of 6
        label_str = f'{LABEL_9CAT.get(cat, cat)} (r={rho:.2f}, {format_p_value(p)})'
        handle = Line2D([], [], marker='o', color='w', markerfacecolor=colors[cat], 
                        markersize=6, markeredgecolor='none', label=label_str)
        handles.append(handle)
        labels.append(label_str)

    # Fitted line
    rho_all, p_all = spearman_stats(x_all, y_all)
    z = np.polyfit(x_all, y_all, 1)
    xs = np.linspace(x_all.min(), x_all.max(), 200)
    ax.plot(xs, np.polyval(z, xs), 'k--', lw=1.4, alpha=0.85, zorder=4)
    
    all_label = f'All (r={rho_all:.2f}, {format_p_value(p_all)})'
    all_handle = Line2D([], [], color='black', linestyle='--', lw=1.4, label=all_label)
    handles.append(all_handle)
    labels.append(all_label)

    ax.set_xlabel('log$_{10}$(genome size / bp)')
    ax.set_ylabel('log$_{10}$(iM count + 1)')
    ax.grid(alpha=0.2, ls='--', lw=0.5); ax.set_axisbelow(True)
    
    ax.legend(handles=handles, labels=labels, loc='upper left', fontsize=6.5, handletextpad=0.3, frameon=False)
    
    if title: ax.set_title(title, fontsize=10, fontweight='bold', loc='left', pad=8)

# D main figure
fig, ax = plt.subplots(figsize=(5.5, 4.5))
draw_D(ax, df, ORDER_4, COLORS_DOMAIN, 'domain4', 'D  Genome size vs iM count')
save(fig, f'{OUT_MAIN}/panel_D_4domain.png')

# D supplementary figure (9 categories)
fig, ax = plt.subplots(figsize=(6.5, 5))
draw_D(ax, df, ORDER_9, COLORS_9CAT, 'Classification_9cat', 'D  Genome size vs iM count (9 cat)')
save(fig, f'{OUT_SUPP}/panel_D_9cat.png')

# D supplementary figure (6 eukaryotic subcategories)
fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
for ax, cat in zip(axes.flat, EUK_6):
    draw_D(ax, df[df['Classification_9cat'] == cat], [cat], COLORS_9CAT, 'Classification_9cat', 
           title=LABEL_9CAT.get(cat, cat)) # titles consistently use LABEL_9CAT
plt.tight_layout()
save(fig, f'{OUT_SUPP}/panel_D_euk6.png')

print('[OK] D saved')
