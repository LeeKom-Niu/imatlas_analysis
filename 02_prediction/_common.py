# _common.py
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

BASE = '/datapool/home/2023200496/niulk/my_project/20_iM_DB_analysis/02_correlation_plots'
DATA_PATH = f'{BASE}/data/iM_summary_clean.tsv'
OUT_MAIN  = f'{BASE}/figures/main_panels'
OUT_SUPP  = f'{BASE}/figures/supplementary'
os.makedirs(OUT_MAIN, exist_ok=True)
os.makedirs(OUT_SUPP, exist_ok=True)

# Fix: viruses (dark yellow) / bacteria (grey-blue) / archaea (grey-green) / eukaryotes (brick red)
# Lowered the saturation to restore the earlier refined look
COLORS_DOMAIN = {
    'Viruses':   '#C29D4A',  # dark gold/earthy yellow
    'Bacteria':  '#3A6B85',  # dark grey-blue
    'Archaea':   '#5C8F50',  # grey-green
    'Eukaryota': '#B65050',  # brick red
}

# 9-category supplementary figure: archaea/bacteria/virus stay consistent; the 6 eukaryotic subcategories use a restrained supporting palette
COLORS_9CAT = {
    'Archaea':          '#5C8F50',  # green
    'Bacteria':         '#3A6B85',  # blue
    'Viral':            '#C29D4A',  # yellow
    'Protozoa':         '#D96C6C',  # soft red
    'Fungi':            '#9467BD',  # purple
    'Plant':            '#8C564B',  # brown
    'Invertebrate':     '#E377C2',  # pink
    'Vertebrate Other': '#7F7F7F',  # grey
    'Mammalian':        '#C0392B',  # dark red
}

# Display labels for the legend and titles
LABEL_9CAT = {
    'Archaea':          'Archaea',
    'Bacteria':         'Bacteria',
    'Viral':            'Viruses',
    'Protozoa':         'Protozoa',
    'Fungi':            'Fungi',
    'Plant':            'Plants',
    'Invertebrate':     'Invertebrates',
    'Vertebrate Other': 'Other vertebrates',
    'Mammalian':        'Mammals',
}

# Internal data mapping
CAT_MAP_4 = {
    'Archaea':'Archaea','Bacteria':'Bacteria','Viral':'Viruses',
    'Protozoa':'Eukaryota','Fungi':'Eukaryota','Plant':'Eukaryota',
    'Invertebrate':'Eukaryota','Vertebrate Other':'Eukaryota','Mammalian':'Eukaryota',
}

ORDER_4 = ['Viruses','Archaea','Bacteria','Eukaryota']
ORDER_9 = ['Archaea','Bacteria','Viral','Protozoa','Fungi','Plant',
           'Invertebrate','Vertebrate Other','Mammalian']
EUK_6 = ['Protozoa', 'Fungi', 'Plant', 'Invertebrate', 'Vertebrate Other', 'Mammalian']

def set_style(fontsize=9):
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.weight': 'normal',
        'font.size': fontsize,
        'axes.labelweight': 'normal',
        'axes.titleweight': 'bold',
        'axes.labelsize': fontsize + 1,
        'axes.titlesize': fontsize + 2,
        'axes.linewidth': 0.9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'xtick.major.width': 0.9,
        'ytick.major.width': 0.9,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'legend.fontsize': fontsize - 1,
        'legend.frameon': False,
        'lines.linewidth': 1.2,
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def load_data():
    df = pd.read_csv(DATA_PATH, sep='\t')
    df['domain4']     = df['Classification_9cat'].map(CAT_MAP_4)
    df['GC_percent']  = df['GC_content'] * 100.0
    df['log10_size']  = np.log10(df['genome_length_bp'].clip(lower=1))
    df['log10_im1']   = np.log10(df['im_count'].clip(lower=0) + 1)
    df['log10_den']   = np.log10(df['density_Mb'].clip(lower=0.01))
    return df

def format_p_value(p):
    if p < 1e-300: return "p < 1e-300"
    if p < 0.0001:
        exp = int(np.floor(np.log10(p)))
        coeff = p / (10 ** exp)
        return f"p = {coeff:.1f}e{exp}"
    elif p < 0.001:
        return f"p < 0.001"
    else:
        return f"p = {p:.3f}"

def spearman_stats(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5: return np.nan, np.nan
    return spearmanr(x[m], y[m])

def marker_size(n):
    if n > 30000: return 1.5
    if n > 5000:  return 3.0
    if n > 500:   return 7
    return 12

def save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches='tight')
    fig.savefig(path.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    plt.close(fig)
