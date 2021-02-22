import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

'''Plotting functions'''

def plot_mt_vs_counts(data, ax, mt_thr=0.5, fontsize=11):    
    # Plot scatter
    ax.scatter(x=data.total_counts, y=data.pct_counts_mt, s=1, c='gray')
    ax.axhline(y=mt_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Fraction MT counts", fontsize=fontsize)


def plot_ngenes_vs_counts(data, ax, gene_thr=6000, fontsize=11):
    # Plot scatter
    ax.scatter(x=data.total_counts, y=data.n_genes_by_counts, s=1, c='gray')
    ax.axhline(y=gene_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Number of genes expr", fontsize=fontsize)


def plot_doublet_scores(data, ax, doublet_thr=0.5, fontsize=11):
    # Plot histogram
    ax.hist(data.doublet_score, bins=100, color='gray')
    ax.axvline(x=doublet_thr, linestyle='--', color="black")
    ax.set_xlabel('Droplet score distribution', fontsize=fontsize)
    

def plot_diss_scores(data, ax, diss_thr=0.5, fontsize=11):
    # Plot histogram
    ax.hist(data.diss_score, bins=100, color='gray')
    ax.axvline(x=diss_thr, linestyle='--', color="black")
    ax.set_xlabel('Dissociation score distribution', fontsize=fontsize)


def plot_ncell_diff(data, ax, labels, n_rem, fontsize=11):
    # Plot Barplot
    for label, n in zip(labels, n_rem):
        ax.bar(label, n)
    ax.set_xlabel('Number of cells removed', fontsize=fontsize)
    ax.tick_params(axis='x', rotation=45)
    
def plot_ngene_diff(adata, ax, fontsize=11):
    ax.set_title('Num genes filtered', fontsize=fontsize)
    ax.bar(x="Before", height=adata.uns['hvg']['ngene'])
    ax.bar(x="After", height=adata.shape[1])
    
def plot_hvg_nbatches(data, ax, fontsize=11):
    for nbatches in np.flip(np.unique(data.highly_variable_nbatches)):
        num = data[data.highly_variable_nbatches == nbatches].shape[0]
        ax.bar(str(nbatches), num, color="gray")
    ax.set_title('Num shared HVG by num samples',fontsize=fontsize)
    
def plot_sorted_rank(data, col, ax, fontsize=11):
    xdim = np.arange(len(data))
    ysort = np.flip(np.sort(data[col]))
    ax.set_title('Ranked {0}'.format(col), fontsize=fontsize)
    ax.plot(xdim, ysort, c='grey')