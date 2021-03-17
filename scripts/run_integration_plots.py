import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from plotting import plot_ngene_diff, plot_hvg_nbatches, plot_ngenes_vs_counts, plot_sorted_rank

'''Plot integrated data set summary plots'''

# Read integrated object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Plot HVG filtering QC plots
fig = plt.figure(figsize=(9,9), dpi=150)
fig.suptitle('HVG filtering QC plots', fontsize=11)
gs = fig.add_gridspec(3, 3)

ax = fig.add_subplot(gs[0,0])
plot_ngene_diff(adata, ax)

ax = fig.add_subplot(gs[0,1])
plot_hvg_nbatches(adata.var, ax)

ax = fig.add_subplot(gs[0,2])
plot_ngenes_vs_counts(adata.obs, ax, gene_thr=np.nan)

ax = fig.add_subplot(gs[1,0:2])
lst_samples = list(adata.obs.sample_id[np.unique(adata.obs.batch.astype(int), return_index=True)[1]])
sc.pl.violin(adata, 'n_genes_by_counts', ax=ax, rotation=45, 
             groupby='sample_id', stripplot=False, show=False, order=lst_samples)

ax = fig.add_subplot(gs[1,2])
plot_sorted_rank(adata.var, 'means', ax)

ax = fig.add_subplot(gs[2,0:2])
sc.pl.violin(adata, 'total_counts', ax=ax, rotation=45, 
             groupby='sample_id', stripplot=False, show=False, order=lst_samples)

ax = fig.add_subplot(gs[2,2])
plot_sorted_rank(adata.var, 'dispersions_norm', ax)

fig.tight_layout()
fig.subplots_adjust(top=0.9)
fig.set_facecolor('white')

# Write to png
fig.savefig('../plots/hvg.png')


# Proj plots
fig, ax = plt.subplots(2,3, figsize=(12,6), dpi=200)
fig.suptitle('Summary projections', fontsize=11)
ax = ax.flatten()
sc.pl.umap(adata, color='sample_id', ax=ax[0], frameon=False, return_fig=False, show=False)
sc.pl.umap(adata, color='condition', ax=ax[1], frameon=False, return_fig=False, show=False)
sc.pl.umap(adata, color='doublet_score', ax=ax[2], frameon=False, return_fig=False, show=False)
sc.pl.umap(adata, color='diss_score', ax=ax[3], frameon=False, return_fig=False, show=False)
sc.pl.umap(adata, color='n_genes_by_counts', ax=ax[4], frameon=False, return_fig=False, show=False)
sc.pl.umap(adata, color='total_counts', ax=ax[5], frameon=False, return_fig=False, show=False)

# Adjust plots
fig.tight_layout()
fig.subplots_adjust(top=0.88)
fig.set_facecolor('white')

# Write to png
fig.savefig('../plots/proj_integration.png')
