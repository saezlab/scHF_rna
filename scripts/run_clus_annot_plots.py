import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

'''Plot the resulting cluster annotations'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)


# Plot top 4 resolution by scores umap projections
fig, axes = plt.subplots(2,2, figsize=(8,6), dpi=300)
axes = axes.flatten()
fig.suptitle('Silhouette scores', fontsize=11)
for i, ax in enumerate(axes):
    res, score = adata.uns['silhouette_scores'][i]
    title = '{0} score:{1:.5}'.format(res, float(score))
    sc.pl.umap(adata, color=res, title=title, frameon=False, ax=ax)

# Adjust plots
fig.tight_layout()