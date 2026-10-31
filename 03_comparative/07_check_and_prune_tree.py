#!/usr/bin/env python3
import pandas as pd
from Bio import Phylo

BASE_DIR = '/datapool/home/2023200496/niulk/my_project/09_type'
TREE_FILE = f'{BASE_DIR}/filtered_tree_unique.nwk'
TABLE_FILE = f'{BASE_DIR}/tree_total_density_filtered.tsv' 
OUT_TREE = f'{BASE_DIR}/final_tree_98sp.nwk'

# 1. Read tree file and extract leaf nodes
print(f"Reading tree file: {TREE_FILE}")
tree = Phylo.read(TREE_FILE, "newick")
tree_leaf_names = set([clade.name for clade in tree.get_terminals() if clade.name])
print(f"Tree file has {len(tree_leaf_names)} leaf nodes.")

# 2. Read table and extract species names
print(f"Reading data table: {TABLE_FILE}")
df = pd.read_csv(TABLE_FILE, sep='\t')
table_names = set(df['Genome'])
print(f"Data table has {len(table_names)} species.")

# 3. Safe tree pruning (traverse all terminal nodes; prune those not in the table)
print("\nStarting safe tree pruning...")
# Must convert to list for traversal
terminals = tree.get_terminals()
for clade in terminals:
    if clade.name not in table_names:
        tree.prune(clade)

# 4. Check the pruned result
final_terminals = tree.get_terminals()
print(f"✅ Pruning complete! Tree now has {len(final_terminals)} species.")

# Save new tree
Phylo.write(tree, OUT_TREE, "newick")
print(f"✅ Perfectly matched tree file saved to:\n{OUT_TREE}")

# 5. Also check whether internal node names carry formats like [100]
# (tvBOT sometimes does not support internal node names with bootstrap values; if the tree still looks odd, run this cleanup)
print("\nCleaning bootstrap values of internal nodes...")
for clade in tree.get_nonterminals():
    if clade.name and '[' in clade.name:
        clade.name = clade.name.split('[')[0] # Strip suffixes like [100]
Phylo.write(tree, f'{BASE_DIR}/final_tree_98sp_clean.nwk', "newick")
print(f"✅ Clean tree file with extraneous internal node text removed saved to:\n{BASE_DIR}/final_tree_98sp_clean.nwk")
