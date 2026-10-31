#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import matplotlib as mpl

# 1. Global style settings (enlarged fonts)
mpl.rcParams.update({
    'font.family': 'DejaVu Sans',  
    'font.size': 11,               # increase base font size
    'axes.titlesize': 15,          # subplot title font
    'axes.labelsize': 14,          # axis label font
    'xtick.labelsize': 12,         # X ticks
    'ytick.labelsize': 12,         # Y ticks
    'axes.linewidth': 0.8,         
    'lines.linewidth': 1.5,        
    'legend.fontsize': 10,         # legend font
    'legend.frameon': False,       
    'pdf.fonttype': 42,            # editable PDF fonts
    'ps.fonttype': 42,
    'svg.fonttype': 'none',        # editable SVG fonts
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

# 2. Paths and colors
BASE_DIR = '/datapool/home/2023200496/niulk/my_project/066_TSS_TES'
OUT_DIR = os.path.join(BASE_DIR, 'Final_Figure5_Plots')
os.makedirs(OUT_DIR, exist_ok=True)

smooth_sigma = 10
window_size = 500
major_xticks = np.arange(-500, 501, 200)

DATA_PATHS = {
    'Archaea':   os.path.join(BASE_DIR, '03_archaea_domain_0506', 'enrichment_results_TSS_TES_1bp_60process_opt'),
    'Bacteria':  os.path.join(BASE_DIR, '04_bacteria_domain_0506', 'enrichment_results_TSS_TES_1bp_60process_opt'),
    'Eukaryota': os.path.join(BASE_DIR, '05_eukaryote_domain_0506', 'enrichment_results_TSS_TES_1bp_eukaryote_60process'),
    'Viruses':   os.path.join(BASE_DIR, '06_viral_domain_0506', 'enrichment_results_TSS_TES_1bp_60process_opt')
}

COLORS_MAIN = {
    'Archaea':   {'template': '#A04000', 'non_template': '#E59866'}, 
    'Bacteria':  {'template': '#154360', 'non_template': '#3498DB'}, 
    'Eukaryota': {'template': '#0E6251', 'non_template': '#2ECC71'}, 
    'Viruses':   {'template': '#4A235A', 'non_template': '#9B59B6'}  
}

COLORS_SUPP = {'template': '#3574B2', 'non_template': '#D43F3A'}

# Change: update names and order
CATEGORY_NAMES = {
    'vertebrate_mammalian': 'Mammals',
    'vertebrate_other': 'Other vertebrates',
    'invertebrate': 'Invertebrates',
    'fungi': 'Fungi',
    'plant': 'Plants',
    'protozoa': 'Protozoa'
}

# 3. Data processing and plotting
def load_data(base_path, region):
    file_path = os.path.join(base_path, f"{region}_combined_results.tsv")
    if not os.path.exists(file_path):
        return None
    
    df = pd.read_csv(file_path, sep="\t")
    
    for col in ['template_enrich', 'non_template_enrich', 'template_lower', 'template_upper', 'non_template_lower', 'non_template_upper']:
        if col in df.columns:
            df[col] = gaussian_filter1d(df[col], sigma=smooth_sigma)
            
    def adjust_ci(df, col_prefix):
        mean = df[f"{col_prefix}_enrich"]
        lower = df[f"{col_prefix}_lower"]
        upper = df[f"{col_prefix}_upper"]
        
        width_half = np.abs(upper - lower) / 2
        new_lower = mean - width_half
        new_upper = mean + width_half
        return new_lower, new_upper

    df["template_lower_adj"], df["template_upper_adj"] = adjust_ci(df, "template")
    df["non_template_lower_adj"], df["non_template_upper_adj"] = adjust_ci(df, "non_template")
    
    df = df.dropna(subset=['template_enrich', 'non_template_enrich'])
    return df

def plot_panel(ax, df, region, colors, title="", letter="", show_ylabel=False, show_legend=False):
    if df is None:
        ax.text(0.5, 0.5, "No Data", ha='center', va='center')
        return

    # Draw confidence intervals
    ax.fill_between(df["position"], df["template_lower_adj"], df["template_upper_adj"], 
                    color=colors['template'], alpha=0.2, edgecolor='none', lw=0)
    ax.fill_between(df["position"], df["non_template_lower_adj"], df["non_template_upper_adj"], 
                    color=colors['non_template'], alpha=0.2, edgecolor='none', lw=0)
    
    # Draw main lines
    ax.plot(df["position"], df["template_enrich"], color=colors['template'], label="Template", lw=1.5)
    ax.plot(df["position"], df["non_template_enrich"], color=colors['non_template'], label="Non-template", lw=1.5)
    
    # Reference lines
    ax.axvline(x=0, color='grey', linestyle='--', lw=0.8)
    ax.axhline(y=1, color='grey', linestyle=':', lw=0.8, alpha=0.8)
    
    # Axis settings (font sizes already enlarged in global rcParams)
    ax.set_xlim(-window_size, window_size)
    ax.set_xticks(major_xticks)
    ax.set_title(title, loc='left', fontweight='bold', pad=5)
    ax.set_xlabel(f'Distance to {region} (bp)')
    
    if show_ylabel:
        ax.set_ylabel('Normalized iM Enrichment')
    
    # Adaptive Y-axis
    y_vals = pd.concat([
        df["template_enrich"], df["non_template_enrich"],
        df["template_lower_adj"], df["template_upper_adj"],
        df["non_template_lower_adj"], df["non_template_upper_adj"]
    ])
    y_min = y_vals.min()
    y_max = y_vals.max()
    
    data_range = max(y_max - y_min, 0.5)
    bottom = max(0, y_min - 0.15 * data_range)
    top = y_max + 0.15 * data_range
    ax.set_ylim(bottom, top)
    
    # Panel letter
    if letter:
        ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontsize=22, fontweight='bold', va='bottom', ha='left')
    
    # Legend
    if show_legend:
        ax.legend(loc='upper right', frameon=False, handlelength=1.2)

# 4. Main program
def main():
    print("Generating main figure (Figure 5 A-H)...")
    fig, axes = plt.subplots(4, 2, figsize=(12, 18)) # slightly larger canvas to prevent crowding with large fonts
    
    domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
    letters = [['A', 'B'], ['C', 'D'], ['E', 'F'], ['G', 'H']]
    
    for i, domain in enumerate(domains):
        df_tss = load_data(DATA_PATHS[domain], 'TSS')
        plot_panel(axes[i, 0], df_tss, 'TSS', COLORS_MAIN[domain], 
                   title=f"{domain} (TSS)", letter=letters[i][0], 
                   show_ylabel=True, show_legend=True)
        
        df_tes = load_data(DATA_PATHS[domain], 'TES')
        plot_panel(axes[i, 1], df_tes, 'TES', COLORS_MAIN[domain], 
                   title=f"{domain} (TES)", letter=letters[i][1], 
                   show_ylabel=False, show_legend=True)
    
    plt.tight_layout()
    main_out = os.path.join(OUT_DIR, 'Figure5_Main_TSS_TES.pdf')
    
    fig.savefig(main_out)
    fig.savefig(main_out.replace('.pdf', '.png'), dpi=300)
    fig.savefig(main_out.replace('.pdf', '.svg'))
    plt.close(fig)
    print(f"✅ Main figure saved: {main_out} (pdf/png/svg)")

    print("\nGenerating eukaryotic subcategory supplementary figures...")
    euk_category_dir = os.path.join(BASE_DIR, '05_eukaryote_domain_0506', 'category_plots_eukaryote')
    
    for cat_key, cat_name in CATEGORY_NAMES.items():
        cat_path = os.path.join(euk_category_dir, cat_key)
        if not os.path.exists(cat_path):
            continue
            
        fig_cat, axes_cat = plt.subplots(1, 2, figsize=(14, 5)) # slightly larger supplementary figure size
        
        df_tss = load_data(cat_path, 'TSS')
        plot_panel(axes_cat[0], df_tss, 'TSS', COLORS_SUPP, 
                   title=f"{cat_name} (TSS)", letter='A', 
                   show_ylabel=True, show_legend=True)
        
        df_tes = load_data(cat_path, 'TES')
        plot_panel(axes_cat[1], df_tes, 'TES', COLORS_SUPP, 
                   title=f"{cat_name} (TES)", letter='B', 
                   show_ylabel=False, show_legend=True)
        
        plt.tight_layout()
        cat_out = os.path.join(OUT_DIR, f'FigS_{cat_key}_TSS_TES.pdf')
        
        fig_cat.savefig(cat_out)
        fig_cat.savefig(cat_out.replace('.pdf', '.png'), dpi=300)
        fig_cat.savefig(cat_out.replace('.pdf', '.svg'))
        plt.close(fig_cat)
        print(f"✅ Subfigure saved: {cat_out} (pdf/png/svg)")

    print("\n🎉 All figures complete! Check the Final_Figure5_Plots folder.")

if __name__ == "__main__":
    main()
