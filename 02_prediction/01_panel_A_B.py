#!/usr/bin/env python3
from _common import *
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

set_style()
df = load_data()

def draw_A(ax, df, order, colors):
    counts = df['Classification_9cat'].value_counts().reindex(
        [o for o in order if o in df['Classification_9cat'].unique()])
    if order == ORDER_4:
        counts = df['domain4'].value_counts().reindex(order)
    cl = [colors[c] for c in counts.index]
    def autopct(p): return f'{p:.1f}%' if p > 4 else ''
    wedges, _, _ = ax.pie(
        counts.values, colors=cl, startangle=90, counterclock=False,
        autopct=autopct, pctdistance=0.80,
        wedgeprops=dict(width=0.42, edgecolor='white', linewidth=1.5),
        textprops=dict(fontsize=8, color='white', fontweight='bold'))
    ax.text(0, 0, f'Total\n{counts.sum():,}', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#333', linespacing=1.4)
    labels = [f'{c} ({counts[c]:,})' for c in counts.index]
    ax.legend(wedges, labels, loc='upper center', bbox_to_anchor=(0.5, -0.02),
              ncol=2, frameon=False, fontsize=7.5)

def draw_B(ax, df, order, colors):
    key = 'domain4' if order == ORDER_4 else 'Classification_9cat'
    sub = df.copy(); sub['grp'] = sub[key]
    stats = sub.groupby('grp').agg(
        n=('Genome','count'),
        total=('im_count','sum'),
        mean_d=('density_Mb','mean'),
        std_d=('density_Mb','std')).reindex(order).dropna(how='all')
    stats['sem'] = stats['std_d'] / np.sqrt(stats['n'])

    colors_list = [colors[c] for c in stats.index]
    x = np.arange(len(stats))
    ax2 = ax.twinx()
    
    # Draw bars and right-axis points
    bars = ax.bar(x, stats['total'], color=colors_list, width=0.6,
                  edgecolor='black', linewidth=0.8, zorder=3)
    ax.set_yscale('log')
    ax.set_ylim(stats['total'].min()*0.3, stats['total'].max()*30)
    ax.minorticks_off()  # Remove minor tick marks from the log axis
    
    for b, v in zip(bars, stats['total']):
        ax.text(b.get_x()+b.get_width()/2, v*1.2, f'{v:,.0f}',
                ha='center', va='bottom', fontsize=7.5)
                
    ax2.errorbar(x, stats['mean_d'], yerr=stats['sem'], fmt='o',
                 markersize=6, capsize=3, color='#222',
                 markerfacecolor='white', markeredgewidth=1.3,
                 elinewidth=1.0, zorder=5)
    ax2.set_ylim(0, stats['mean_d'].max()*1.5)
    ax2.minorticks_off()
    
    for xi, m, s in zip(x, stats['mean_d'], stats['sem']):
        ax2.annotate(f'{m:.1f}', (xi, m+s), xytext=(0, 7),
                     textcoords='offset points', ha='center', fontsize=7.5)

    # Show the right axis and thicken the spine
    ax2.spines['right'].set_visible(True)
    ax2.spines['right'].set_linewidth(0.9)
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(axis='y', width=0.9, length=3.5)

    ax.set_xticks(x)
    ax.set_xticklabels(stats.index, fontsize=8.5, rotation=30 if len(stats) > 6 else 0,
                       ha='right' if len(stats) > 6 else 'center')
    ax.set_ylabel('Total iM count (log)', fontsize=9.5)
    ax2.set_ylabel('Mean density (iMs/Mb)', fontsize=9.5)
    ax.grid(axis='y', ls='--', lw=0.5, alpha=0.25); ax.set_axisbelow(True)
    
    ax.legend(handles=[Patch(facecolor='#888', edgecolor='k', label='Total iMs (log)'),
               Line2D([0],[0],marker='o',color='#222',markerfacecolor='white',
                      markersize=6, markeredgewidth=1.3, lw=0,
                      label='Mean density ± SEM')],
              loc='upper left', fontsize=7)

for tag, order, colors in [('4domain', ORDER_4, COLORS_DOMAIN), ('9cat', ORDER_9, COLORS_9CAT)]:
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    draw_A(ax, df, order, colors)
    save(fig, f'{OUT_MAIN if tag=="4domain" else OUT_SUPP}/panel_A_{tag}.png')

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    draw_B(ax, df, order, colors)
    save(fig, f'{OUT_MAIN if tag=="4domain" else OUT_SUPP}/panel_B_{tag}.png')

fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 4.5), gridspec_kw=dict(wspace=0.35, width_ratios=[1, 1.15]))
draw_A(axL, df, ORDER_4, COLORS_DOMAIN)
draw_B(axR, df, ORDER_4, COLORS_DOMAIN)
save(fig, f'{OUT_MAIN}/panel_AB_combined.png')
print('[OK] A, B saved')
