import scanpy as sc
import numpy as np
import pandas as pd

import dorothea
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

# Get dorothea model
model = dorothea.load_regulons()

# Compute activities
dorothea.run(adata, model, center=True, norm=True, scale=True, num_perm=100, use_raw=False)

dfs = []
for cell_type in cell_types:
    # Filter by cell type
    subadata = adata[(adata.obs['cell_type'] == cell_type)]
    
    # Compute DAP
    subadata.obs['condition'] = subadata.obs['condition'].astype(str)
    design = get_design(subadata.obs, 'condition')
    contr_matrix = get_contrast(design, contr_dict)
    data = dorothea.extract(subadata)
    data = pd.DataFrame(data.X.T, index=data.var.index, columns=data.obs.index)
    df = limma_fit(data, design, contr_matrix).sort_values(['contrast', 'pvals'])
    df['cell_type'] = cell_type
    dfs.append(df)

# Merge all dfs
df = pd.concat(dfs)

# Save
os.makedirs('../plot_data/func/', exist_ok=True)
df.to_csv('../plot_data/func/dorothea.csv', index=False)
adata.obsm['dorothea'].to_csv('../plot_data/func/dorothea_act.csv', index=True, header=True)