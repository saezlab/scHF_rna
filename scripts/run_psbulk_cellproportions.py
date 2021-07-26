import pandas as pd
import numpy as np

from scipy import stats
from statsmodels.stats.multitest import multipletests

# Read cell proportions
df = pd.read_csv('../plot_data/var_shifts/cell_proportions.csv')

# Iterate per cell type and condition
table = []
for cell_type in np.unique(df['CellType']):
    for condition in np.unique(df['Condition']):
        group_a = df[(df['CellType'] == cell_type) & (df['Condition'] == condition)]['Prop'].values
        group_b = df[(df['CellType'] == cell_type) & (df['Condition'] != condition)]['Prop'].values
        
        # Wilcoxon rank-sum statistic
        stat, pval = stats.ranksums(group_a, group_b)
        
        table.append([cell_type, condition, stat, pval])
        
table = pd.DataFrame(table, columns=['cell_type', 'condition', 'stat', 'pval'])

# Correct by FDR
_, pvals_adj, _, _ = multipletests(table['pval'].values, alpha=0.05, method='fdr_bh')
table['adj_pval'] = pvals_adj
table = table.sort_values('pval')

table.to_csv('../tables/cellproportions.csv', index=False)