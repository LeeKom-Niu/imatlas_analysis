#!/usr/bin/env python3
import pandas as pd
from Bio import Phylo

BASE_DIR = '/datapool/home/2023200496/niulk/my_project/09_type'
TREE_FILE = f'{BASE_DIR}/filtered_tree_unique.nwk'
TABLE_FILE = f'{BASE_DIR}/tree_total_density_filtered.tsv' 
OUT_TREE = f'{BASE_DIR}/final_tree_98sp.nwk'

# 1. Read tree file and extract leaf nodes
print(f"正在读取树文件: {TREE_FILE}")
tree = Phylo.read(TREE_FILE, "newick")
tree_leaf_names = set([clade.name for clade in tree.get_terminals() if clade.name])
print(f"树文件中共有 {len(tree_leaf_names)} 个叶子节点。")

# 2. Read table and extract species names
print(f"正在读取数据表: {TABLE_FILE}")
df = pd.read_csv(TABLE_FILE, sep='\t')
table_names = set(df['Genome'])
print(f"数据表中共有 {len(table_names)} 个物种。")

# 3. Safe tree pruning (traverse all terminal nodes; prune those not in the table)
print("\n开始安全修剪树文件...")
# Must convert to list for traversal
terminals = tree.get_terminals()
for clade in terminals:
    if clade.name not in table_names:
        tree.prune(clade)

# 4. Check the pruned result
final_terminals = tree.get_terminals()
print(f"✅ 修剪完成！现在树上有 {len(final_terminals)} 个物种。")

# Save new tree
Phylo.write(tree, OUT_TREE, "newick")
print(f"✅ 完美匹配的树文件已保存至:\n{OUT_TREE}")

# 5. Also check whether internal node names carry formats like [100]
# (tvBOT sometimes does not support internal node names with bootstrap values; if the tree still looks odd, run this cleanup)
print("\n正在清洗内部节点的 bootstrap 值...")
for clade in tree.get_nonterminals():
    if clade.name and '[' in clade.name:
        clade.name = clade.name.split('[')[0] # Strip suffixes like [100]
Phylo.write(tree, f'{BASE_DIR}/final_tree_98sp_clean.nwk', "newick")
print(f"✅ 去除内部节点杂字的纯净版树文件已保存至:\n{BASE_DIR}/final_tree_98sp_clean.nwk")
