import scanpy as sc
import numpy as np
import pandas as pd

import os

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path) 

# Get unique conditions and remove healthy
conditions = np.unique(adata.obs['condition'])

# Get unique cell types
cell_types = np.unique(adata.obs['cell_type'])

# For each condition and cell type compute DEG and store them in a df
dfs = []
for i,cond_a in enumerate(conditions):
    for j,cond_b in enumerate(conditions):
        if j < i:
            for cell_type in cell_types:
                # Filter by cell type
                subadata = adata[(adata.obs['cell_type'] == cell_type)]

                # Drop genes that have 0 expression
                subadata = subadata[:,~(subadata.X == 0).any(axis=0)]

                # Skip if only one sample or less are available
                if np.sum(subadata.obs['condition'] == cond_a) <= 1 or \
                   np.sum(subadata.obs['condition'] == cond_b) <= 1:
                    continue

                # Compute DEG
                sc.tl.rank_genes_groups(subadata, 'condition', groups=[cond_a, cond_b],
                                        reference=cond_a, method='t-test')

                # Get df and store it
                df = sc.get.rank_genes_groups_df(subadata, cond_b)
                df['cell_type'] = cell_type
                df['contrast'] = '{0}-{1}'.format(cond_b, cond_a)
                df.drop('scores', axis=1, inplace=True)
                dfs.append(df)

# Merge all dfs
df = pd.concat(dfs)

# Save
os.makedirs('../plot_data/deg/', exist_ok=True)
df.to_csv('../plot_data/deg/psbulk_deg.csv', index=False)