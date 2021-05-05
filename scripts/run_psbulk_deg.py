import scanpy as sc
import numpy as np
import pandas as pd

import os

from utils import get_design, get_contrast, limma_fit

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path) 

# Get unique conditions
conditions = np.unique(adata.obs['condition'])

# Get unique cell types
cell_types = np.unique(adata.obs['cell_type'])

### Contrasts to test ###

contr_dict = {
    'hf-healthy' : ['hf','healthy'],
    'hf_ckd-healthy' : ['hf_ckd','healthy'],
    'hf_ckd-hf' : ['hf_ckd','hf']
}

###

dfs = []
for cell_type in cell_types:
    # Filter by cell type
    subadata = adata[(adata.obs['cell_type'] == cell_type)]

    # Drop genes that have 0 expression
    subadata = subadata[:,~(subadata.X == 0).any(axis=0)]

    # Compute DEG
    subadata.obs['condition'] = subadata.obs['condition'].astype(str)
    design = get_design(subadata.obs, 'condition')
    contr_matrix = get_contrast(design, contr_dict)
    data = pd.DataFrame(subadata.X.T, index=subadata.var.index, columns=subadata.obs.index)
    deg = limma_fit(data, design, contr_matrix).sort_values(['contrast', 'pvals'])

    deg['cell_type'] = cell_type
    dfs.append(deg)


# Merge all dfs
df = pd.concat(dfs)

# Save
os.makedirs('../plot_data/deg/', exist_ok=True)
df.to_csv('../plot_data/deg/psbulk_deg.csv', index=False)