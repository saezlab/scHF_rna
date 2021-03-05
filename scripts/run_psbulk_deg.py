import scanpy as sc
import numpy as np
import pandas as pd

import os

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path) 

# Get conditions and remove healthy
conditions = adata.obs['condition'][adata.obs['condition'] != 'healthy']

# For each condition and cell type compute DEG and store them in a df
dfs = []
for condition in np.unique(conditions):
    for ctype in np.unique(adata.obs['cell_type']):
        # Filter by cell type
        subadata = adata[(adata.obs['cell_type'] == ctype)]
        # Drop conditions with only one sample
        subadata = subadata[subadata.obs["condition"].duplicated(keep=False)]
        # Drop genes that have 0 expression
        subadata = subadata[:,~(subadata.X == 0).any(axis=0)]

        # Skip if no samples are available
        if np.sum(subadata.obs['condition'] == condition) == 0:
            continue

        # Compute DEG
        sc.tl.rank_genes_groups(subadata, 'condition', reference='healthy', method='t-test')

        # Get df and store it
        df = sc.get.rank_genes_groups_df(subadata, condition)
        df['cell_type'] = ctype
        df['condition'] = condition
        df.drop('scores', axis=1, inplace=True)
        dfs.append(df)

# Merge all dfs
df = pd.concat(dfs)

# Save
os.makedirs('../plot_data/deg/', exist_ok=True)
df.to_csv('../plot_data/deg/deg.csv', index=False)