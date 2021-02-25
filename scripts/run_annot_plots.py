import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

'''Plot the resulting cluster annotations'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Read marker genes
markers = pd.read_csv('../data/markers.csv')
msk = np.array([marker in adata.var.index for marker in markers.Gene])
markers = markers[msk]
marker_genes_dict = dict()
for ctype in np.unique(markers.Type):
    genes = list(markers.Gene[markers.Type == ctype])
    marker_genes_dict[ctype] = genes

# Plot cell type annotation
fig = plt.figure(figsize=(12,12), dpi=150)
fig.suptitle('Cell annotation', fontsize=15)
gs = fig.add_gridspec(3,2)

ax = fig.add_subplot(gs[0:2,:])
sc.pl.umap(adata, color='cell_type', ax=ax, show=False, 
           return_fig=False, frameon=False, legend_loc='on data', size=10, title="")

ax = fig.add_subplot(gs[2,:])
sc.pl.dotplot(adata, marker_genes_dict, 'cell_type', dendrogram=True, ax=ax, show=False, return_fig=False)

fig.tight_layout()
fig.subplots_adjust(top=0.95)
fig.set_facecolor('white')

# Save
fig.savefig('../plots/proj_annotation.png')