#!/usr/bin/env python3
from _common import *
import matplotlib.image as mpimg

# Load each panel image and stitch them into a 4x2 layout
panels = [
    ('A', f'{OUT_MAIN}/panel_A_4domain.png'),
    ('B', f'{OUT_MAIN}/panel_B_4domain.png'),
    ('C', f'{OUT_MAIN}/panel_C_4domain.png'),
    ('D', f'{OUT_MAIN}/panel_D_4domain.png'),
    ('E', f'{OUT_MAIN}/panel_E.png'),
    ('F', f'{OUT_MAIN}/panel_F.png'),
    ('G', f'{OUT_MAIN}/panel_G.png'),
    ('H', f'{OUT_MAIN}/panel_H.png'),
]

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
for ax, (tag, path) in zip(axes.flat, panels):
    img = mpimg.imread(path)
    ax.imshow(img)
    ax.axis('off')
    ax.text(-0.05, 1.02, tag, transform=ax.transAxes,
            fontsize=18, fontweight='bold', va='bottom', ha='left')

plt.tight_layout()
fig.savefig(f'{OUT_MAIN}/figure_main_A_to_H.png',
            dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT_MAIN}/figure_main_A_to_H.pdf',
            dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'[OK] combined main figure -> {OUT_MAIN}/figure_main_A_to_H.png')
