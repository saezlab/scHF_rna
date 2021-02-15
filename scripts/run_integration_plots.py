import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

'''Plot integrated data set summary plots'''

# Read integrated object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Summary plots
fig, ax = plt.subplots(2,2, figsize=(8,6), dpi=300)
fig.suptitle('Summary projections', fontsize=11)
ax = ax.flatten()
sc.pl.umap(adata, color='sample_id', ax=ax[0], frameon=False)
sc.pl.umap(adata, color='condition', ax=ax[1], frameon=False)
sc.pl.umap(adata, color='doublet_score', ax=ax[2], frameon=False)
sc.pl.umap(adata, color='pct_counts_mt', ax=ax[3], frameon=False)

# Adjust plots
fig.tight_layout()
fig.set_facecolor('white')

# Write to png
fig.savefig('../plots/proj_integration.png')
