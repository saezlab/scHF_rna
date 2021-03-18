import scanpy as sc
import numpy as np
import pandas as pd

import dorothea
import os
from utils import lm, rank_func_feature

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

cell_types = np.unique(adata.obs['cell_type'])
conditions = np.unique(adata.obs['condition'])


func = 'dorothea'
# Read Dorothea Regulons for Human
dorothea_hs = dorothea.load_regulons(['A', 'B', 'C'])

results = []
for cell_type in cell_types:
    subadata = adata[adata.obs['cell_type'] == cell_type]
    dorothea.run_scira(subadata, dorothea_hs, norm='c', scale=True)
    names = subadata.obsm[func].columns
    for i,cond_a in enumerate(conditions):
        for j,cond_b in enumerate(conditions):
            if cond_a != cond_b:
                for name in names:
                    coeff, pval = rank_func_feature(subadata, name, cond_a, cond_b, func)
                    result = [cell_type, cond_a, cond_b, name, coeff, pval]
                    results.append(result)

# Transform to df
results = pd.DataFrame(results, columns=['cell_type', 'ref', 'cond', 'name', 'coeff', 'pvals']).sort_values('pvals')

# Save
os.makedirs('../plot_data/func/', exist_ok=True)
results.to_csv('../plot_data/func/dorothea.csv', index=False)