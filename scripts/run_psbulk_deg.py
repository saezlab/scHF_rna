import scanpy as sc
import numpy as np
import pandas as pd

import os

from utils import edgeR
from statsmodels.stats.multitest import multipletests
from scipy.stats import f_oneway

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path) 

# Get unique conditions
conditions = np.unique(adata.obs['condition'])

# Get unique cell types
cell_types = np.unique(adata.obs['cell_type'])

### Contrasts to test ###

contr_dict = {
    'HF-AvsHealthy' : ['HF-A','Healthy'],
    'HF-CKDvsHealthy' : ['HF-CKD','Healthy'],
    'HF-CKDvsHF-A' : ['HF-CKD','HF-A']
}

###

dfs = []
for cell_type in cell_types:
    # Filter by cell type
    subadata = adata[(adata.obs['cell_type'] == cell_type)]
    msk = np.sum(subadata.X, axis=0) > 0
    subadata = subadata[:, msk]
    subadata.obs['condition'] = subadata.obs['condition'].astype(str)
    
    # Subset genes based on ANOVA
    num_g = subadata.shape[1]
    g_msk = np.full(shape=(num_g,), fill_value=False)
    for i in range(num_g):
        X = subadata[:,i].X.flatten()
        groups = [X[subadata.obs['condition'] == condition] for condition in conditions]
        p_value = f_oneway(*groups).pvalue
        if p_value < 0.05:
            g_msk[i] = True
    subadata = subadata.raw.to_adata()[:,msk][:,g_msk]
    print(cell_type, subadata.shape)
    
    # Format data for edgeR
    data = pd.DataFrame(subadata.X.T, index=subadata.var.index, columns=subadata.obs.index)
    meta = subadata.obs
    
    # Compute DEG
    for contrast in contr_dict.keys():
        cond_a, cond_b = contrast.split('vs')
        msk_a = meta['condition']==cond_a
        msk_b = meta['condition']==cond_b
        # Skip if groups have less than 2 samples each
        if np.sum(msk_a) < 2 or np.sum(msk_b) < 2:
            continue
        msk = msk_a | msk_b
        deg = edgeR(data.loc[:,msk], meta.loc[msk], contrast)
        deg['contrast'] = contrast
        deg['cell_type'] = cell_type
        dfs.append(deg)

# Merge all dfs
df = pd.concat(dfs)

# Save
os.makedirs('../plot_data/deg/', exist_ok=True)
df.to_csv('../plot_data/deg/psbulk_deg.csv', index=False)