import numpy as np
import pandas as pd
import pickle

import os
from plotting import *

"""
Script to plot different QC metrics after filtering the data.
"""

# Load meta data
meta = pd.read_csv('../data/metadata.csv')
samples = np.unique(meta['sample_id'])

# Create plots dir
os.makedirs('../plots/', exist_ok=True)

# Initialize sumary variables
total_df = pd.DataFrame(columns=['n_genes_by_counts','total_counts',
                                   'pct_counts_mt','doublet_score'])
total_subsets = np.zeros(7).astype(np.int)
summary_df = []
total_n_rem = np.zeros(5)

# Run QC plots for each sample and store summary
for sample in samples:
    print(sample)
    input_path = '../plot_data/qc/{0}.pkl'.format(sample)
    plot_data = pickle.load(open(input_path, "rb"))

    # Filter params
    mt_thr = plot_data['mt_thr']
    gene_qnt = plot_data['gene_qnt']
    doublet_thr = plot_data['doublet_thr']
    diss_qnt = plot_data['diss_qnt']
    df = plot_data['df']
    gene_thr = np.quantile(df.n_genes_by_counts, gene_qnt)
    diss_thr = np.quantile(df.diss_score, diss_qnt)

    # Create QC plots
    fig, axes = plt.subplots(2,3, figsize=(11,8.5), dpi=150)
    fig.suptitle('QC metrics {0}'.format(sample), fontsize=11)
    axes = axes.flatten()

    # Plot MT
    plot_mt_vs_counts(df, axes[0], mt_thr=mt_thr)

    # Plot ngenes
    plot_ngenes_vs_counts(df, axes[1], gene_thr=gene_thr)

    # Plot lost cells
    labels = ['MT', 'Gene', 'Doublet', 'Diss', 'Total']
    msks = np.array([
        df.pct_counts_mt > mt_thr,
        df.n_genes_by_counts > gene_thr,
        df.doublet_score > doublet_thr,
        df.diss_score > diss_thr
    ])
    msks = np.vstack((msks, [np.sum(msks, axis=0) > 0]))
    n_rem = np.sum(msks, axis=1)
    total_n_rem += n_rem
    plot_ncell_diff(df, axes[2], labels=labels, n_rem=n_rem)
    
    # Plot doublet score
    plot_doublet_scores(df, axes[3], doublet_thr=doublet_thr)

    # Plot diss score
    plot_diss_scores(df, axes[4], diss_thr=diss_thr)
    
    # Filter
    msk = (df.n_genes_by_counts < gene_thr) & \
      (df.pct_counts_mt < mt_thr) &  \
      (df.doublet_score < doublet_thr) & \
      (df.diss_score < diss_thr)
    df = df[msk]

    # Plot diff cells
    pd.DataFrame([[sample, len(msk), sum(msk)]], 
                 columns=["Sample ID", "Before", "After"]).plot.bar(x=0, width=0.5, ax=axes[5])
    axes[5].set_title('Lost cells', fontsize=11)
    
    # Adjust plots
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    fig.set_facecolor('white')

    # Write to png
    fig.savefig('../plots/qc_{0}.png'.format(sample))

    # Append
    total_df = total_df.append(df, ignore_index=True)
    summary_df.append([sample, len(msk), sum(msk)])
        
# Summary plots
fig, axes = plt.subplots(2,3, figsize=(11,8.5), dpi=150)
fig.suptitle('Merged QC metrics', fontsize=11)
axes = axes.flatten()

plot_mt_vs_counts(total_df, axes[0], mt_thr=np.nan)

plot_ngenes_vs_counts(total_df, axes[1], gene_thr=np.nan)

plot_ncell_diff(df, axes[2], labels=labels, n_rem=total_n_rem)

plot_doublet_scores(total_df, axes[3], doublet_thr=np.nan)

plot_diss_scores(total_df, axes[4], diss_thr=np.nan)

pd.DataFrame(summary_df, columns=["Sample ID", "Before", "After"]).plot.bar(x=0, ax=axes[5])
axes[5].set_title('Lost cells', fontsize=11)

fig.tight_layout()
fig.subplots_adjust(top=0.88)
fig.set_facecolor('white')

# Save
fig.savefig('../plots/qc_summary.png')
        