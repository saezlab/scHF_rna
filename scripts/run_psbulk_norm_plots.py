import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

# Open AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

# Init figure
fig, axes = plt.subplots(3, 4, figsize=(4*4, 3*4), dpi=150, sharey=True, sharex=True)
axes = axes.flatten()
ax = axes[0]

# Set 0 to nan to not bias violins
adata.X[adata.X == 0] = np.nan

# Plot norm for all cell types
v_x = np.repeat(adata.obs['sample_id'].tolist(), adata.shape[1])
v_y = adata.X.flatten()
v_hue = np.repeat(adata.obs['cell_type'].tolist(), adata.shape[1])

ax.set_title('All')
sns.violinplot(x=v_x, y=v_y, ax=ax, inner="quartile")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylabel('Norm expr')

# Plot norm for specific cell types
for ctype,ax in zip(np.unique(adata.obs['cell_type']), axes[1:]):
    subadata = adata[adata.obs['cell_type'] == ctype].copy()

    v_x = np.repeat(subadata.obs['sample_id'].tolist(), subadata.shape[1])
    v_y = subadata.X.flatten()
    v_hue = np.repeat(subadata.obs['cell_type'].tolist(), subadata.shape[1])
    
    ax.set_title(ctype)
    sns.violinplot(x=v_x, y=v_y, ax=ax, inner="quartile")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    
fig.tight_layout()
fig.set_facecolor('white')

fig.savefig('../plots/psbulk_qc.png')