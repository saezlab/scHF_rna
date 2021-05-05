import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

import scanpy as sc

from plotting import violins 

# Read data
input_path = '../qc_data/bulk.h5ad'
adata = sc.read_h5ad(input_path).raw.to_adata()
adata.obs['sample_id'] = adata.obs.index.values
adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
data = pd.DataFrame(adata.X.T, index=adata.var.index, columns=adata.obs.index)
data.columns.name = None

# Plot
fig = plt.figure(figsize=(9,9), tight_layout=True, dpi=150, facecolor='white')
fig.suptitle('Bulk QC plots', fontsize=11)
gs = fig.add_gridspec(2, 2)

ax = fig.add_subplot(gs[0,:])
violins(data, adata.obs, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

ax = fig.add_subplot(gs[1,0])
sc.pl.pca(adata, color='condition', ax=ax, return_fig=False, show=False, size=200)

ax = fig.add_subplot(gs[1,1])
sc.pl.pca(adata, color='sample_id', ax=ax, return_fig=False, show=False, size=0, legend_loc='on data', legend_fontweight='light')

# Save
fig.savefig('../plots/bulk_qc.png')