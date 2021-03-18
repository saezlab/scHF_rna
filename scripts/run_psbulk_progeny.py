import scanpy as sc
import numpy as np
import pandas as pd

import progeny
import os
from utils import lm, rank_func_feature

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

cell_types = np.unique(adata.obs['cell_type'])
conditions = np.unique(adata.obs['condition'])

# Define functional name progeny
func = 'progeny'
results = []
for cell_type in cell_types:
    # SUbset by cell type and run progeny
    subadata = adata[adata.obs['cell_type'] == cell_type]
    progeny.run(subadata, scale=True, organism='Human', top=100)
    names = subadata.obsm[func].columns
    for i,cond_a in enumerate(conditions):
        for j,cond_b in enumerate(conditions):
            if cond_a != cond_b:
                for name in names:
                    # For each pathway, compute significant differences in mean (linear model)
                    coeff, pval = rank_func_feature(subadata, name, cond_a, cond_b, func)
                    result = [cell_type, cond_a, cond_b, name, coeff, pval]
                    results.append(result)

# Transform to df
results = pd.DataFrame(results, columns=['cell_type', 'ref', 'cond', 'name', 'coeff', 'pvals']).sort_values('pvals')

# Save
os.makedirs('../plot_data/func/', exist_ok=True)
results.to_csv('../plot_data/func/progeny.csv', index=False)