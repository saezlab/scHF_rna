import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

'''Plot the resulting cluster annotations'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Read markers genes
markers = pd.read_csv('../data/markers.csv')
msk = np.array([marker in adata.var.index for marker in markers.Gene])
markers = markers[msk]
marker_genes_dict = dict()
for ctype in np.unique(markers.Type):
    genes = list(markers.Gene[markers.Type == ctype])
    marker_genes_dict[ctype] = genes

# List of resolutions
ress = [name for name in adata.obs.columns if name.startswith('leiden_res_')]

for res in ress:
    # Plot projection and dotplot
    fig = plt.figure(figsize=(12,12), dpi=150)
    gs = fig.add_gridspec(3,2)

    ax = fig.add_subplot(gs[0:2,:])
    sc.pl.umap(adata, color=res, ax=ax, show=False, 
               return_fig=False, frameon=False, legend_loc='on data', size=10)

    ax = fig.add_subplot(gs[2,:])
    sc.pl.dotplot(adata, marker_genes_dict, res, dendrogram=True, ax=ax, show=False, return_fig=False)

    fig.tight_layout()
    fig.set_facecolor('white')
    
    # Save
    fig.savefig('../plots/{0}.png'.format(res))
    
