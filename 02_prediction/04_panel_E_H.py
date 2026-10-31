#!/usr/bin/env python3
from _common import *
from scipy.stats import linregress
import matplotlib.image as mpimg

set_style()
df = load_data()

def draw_and_save_EH(df, cat, color, tag, title):
    fig, ax = plt.subplots(figsize=(6, 5))
    sub = df[(df['domain4'] == cat) & (df['density_Mb'] > 0)].copy() if 'domain4' in df else df[df['Classification_9cat'] == cat]
    
    if len(sub) < 5: return
    x = sub['log10_size'].values; y = np.log10(sub['density_Mb']).values
    rho, p = spearman_stats(x, y)
    slope, intercept, _, _, _ = linregress(x, y)

    ax.scatter(x, y, color=color, s=20, alpha=0.6, edgecolor='none', label=title)
    x_min, x_max = x.min() - 0.5, x.max() + 0.5
    ax.plot(np.linspace(x_min,x_max,100), slope*np.linspace(x_min,x_max,100)+intercept, 
            'k--', lw=1.5, label=f'r = {rho:.2f}, {format_p_value(p)}')
    
    ax.set_xlabel('log$_{10}$(Genome size)')
    ax.set_ylabel('log$_{10}$(iM density)')
    ax.legend(loc='upper left')
    ax.set_xlim(x_min, x_max); ax.set_ylim(y.min() - 0.3, y.max() + 0.3)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_title(f'{tag}  {title}', loc='left')
    
    # Save the individual panel
    save(fig, f'{OUT_MAIN}/panel_{tag}.png')
    print(f'[OK] saved {tag} {cat}')

# 1. Generate the four individual panels (E, F, G, H)
# Note: the key names have been updated to the correct Eukaryota and Viruses
cats = [('Archaea', 'E'), ('Bacteria', 'F'), ('Eukaryota', 'G'), ('Viruses', 'H')]
for cat, tag in cats:
    display_title = LABEL_9CAT.get(cat, cat)
    draw_and_save_EH(df, cat, COLORS_DOMAIN[cat], tag, display_title)

# 2. Combine the four E-H panels into a 2x2 composite
print("Combining E-H panels into 2x2...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
# The panels here are used to read the generated images; their order must match the file names
panels = [('E', 'Archaea'), ('F', 'Bacteria'), ('G', 'Eukaryota'), ('H', 'Viruses')]

for ax, (tag, cat) in zip(axes.flat, panels):
    img_path = f'{OUT_MAIN}/panel_{tag}.png'
    try:
        img = mpimg.imread(img_path)
        ax.imshow(img)
        ax.axis('off')
        ax.text(-0.05, 1.05, tag, transform=ax.transAxes,
                fontsize=24, fontweight='bold', va='bottom', ha='left')
    except FileNotFoundError:
        ax.text(0.5, 0.5, f'Missing {tag}', ha='center', va='center')
        ax.axis('off')

plt.tight_layout()
save(fig, f'{OUT_MAIN}/panel_E_to_H_combined_2x2.png')
print(f'[OK] Combined 2x2 saved to {OUT_MAIN}/panel_E_to_H_combined_2x2.png')

# 3. Combine the four E-H panels into a 1x4 composite
print("Combining E-H panels into 1x4...")
fig, axes = plt.subplots(1, 4, figsize=(24, 5)) 

for ax, (tag, cat) in zip(axes.flat, panels):
    img_path = f'{OUT_MAIN}/panel_{tag}.png'
    try:
        img = mpimg.imread(img_path)
        ax.imshow(img)
        ax.axis('off')
        ax.text(-0.05, 1.05, tag, transform=ax.transAxes,
                fontsize=24, fontweight='bold', va='bottom', ha='left')
    except FileNotFoundError:
        ax.text(0.5, 0.5, f'Missing {tag}', ha='center', va='center')
        ax.axis('off')

plt.tight_layout()
save(fig, f'{OUT_MAIN}/panel_E_to_H_combined_1x4.png')
print(f'[OK] Combined 1x4 saved to {OUT_MAIN}/panel_E_to_H_combined_1x4.png')

# 4. Generate the E-H supplementary figure (6 eukaryotic subcategories)
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
for ax, cat in zip(axes.flat, EUK_6):
    sub = df[(df['Classification_9cat'] == cat) & (df['density_Mb'] > 0)].copy()
    if len(sub) > 5:
        x = sub['log10_size'].values; y = np.log10(sub['density_Mb']).values
        rho, p = spearman_stats(x, y)
        slope, intercept, _, _, _ = linregress(x, y)
        ax.scatter(x, y, color=COLORS_9CAT[cat], s=15, alpha=0.6, edgecolor='none')
        x_min, x_max = x.min() - 0.5, x.max() + 0.5
        ax.plot(np.linspace(x_min,x_max,100), slope*np.linspace(x_min,x_max,100)+intercept, 
                'k--', lw=1.5, label=f'r = {rho:.2f}, {format_p_value(p)}')
        
        # Change: use the LABEL_9CAT mapping for titles to show the correct plural names
        ax.set_title(LABEL_9CAT.get(cat, cat)) 
        
        ax.set_xlabel('log$_{10}$(Genome size)')
        ax.set_ylabel('log$_{10}$(iM density)')
        ax.set_xlim(x_min, x_max)
        ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
save(fig, f'{OUT_SUPP}/panel_EH_euk6.png')
print('[OK] saved E-H euk6')
