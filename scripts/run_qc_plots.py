import numpy as np
import pandas as pd
import pickle

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn3_unweighted

"""
Script to plot different QC metrics after filtering the data.
"""


def plot_mt_vs_counts(data, ax, mt_thr=0.5, fontsize=11):    
    # Plot scatter
    ax.scatter(x=data.total_counts, y=data.pct_counts_mt, s=1, c='gray')
    ax.axhline(y=mt_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Fraction mitocondrial counts", fontsize=fontsize)


def plot_ngenes_vs_counts(data, ax, gene_thr=6000, fontsize=11):
    # Plot scatter
    ax.scatter(x=data.total_counts, y=data.n_genes_by_counts, s=1, c='gray')
    ax.axhline(y=gene_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Number of genes expressed", fontsize=fontsize)


def plot_doublet_scores(data, ax, doublet_thr=0.5, fontsize=11):
    # Plot histogram
    ax.hist(data.doublet_score, bins=100, color='gray')
    ax.axvline(x=doublet_thr, linestyle='--', color="black")
    ax.set_xlabel('Droplet score distribution', fontsize=fontsize)
    ax.set_xlim(-0.05,1)
    ax.set_xticks(np.arange(0,1,0.1))


def plot_ncell_diff(data, ax, mt_thr=0.5, gene_thr=6000, doublet_thr=0.5, fontsize=11):
    msk = (df.n_genes_by_counts < gene_thr) & \
          (df.pct_counts_mt < mt_thr) &  \
          (df.doublet_score < doublet_thr)
    # Plot Barplot
    ax.bar(x=["Before","After"],
        height=[data.shape[0], data[msk].shape[0]],
        color="gray"
       )
    ax.set_ylabel("Number of cells", fontsize=fontsize)


def compute_overlaps(data, mt_thr=0.5, gene_thr=6000, doublet_thr=0.5):
    # Build masks
    msk_gene = data.n_genes_by_counts > gene_thr
    msk_mt = data.pct_counts_mt > mt_thr
    msk_doublet = data.doublet_score > doublet_thr
    
    # Count overlaps
    Abc = sum(msk_gene & ~(msk_mt | msk_mt))
    aBc = sum(msk_mt & ~(msk_gene | msk_doublet))
    ABc = sum(msk_gene & msk_mt & ~(msk_doublet))
    abC = sum(msk_doublet & ~(msk_gene | msk_mt))
    AbC = sum(msk_gene & msk_doublet & ~(msk_mt))
    aBC = sum(msk_mt & msk_doublet & ~(msk_gene))
    ABC = sum(msk_gene & msk_mt & msk_doublet)

    return [Abc, aBc, ABc, abC, AbC, aBC, ABC]
    
    
def plot_venn(subsets, ax):
    # Plot Venn
    venn3_unweighted(subsets = subsets, 
          set_labels = ('N. Genes', 'MT', 'Doublet'),
          ax=ax
         )


meta = {
    'healthy' : ["CK114","CK115","CK139","CK140"],
    'acidosis' : ["CK128"],
    'hf' : ["CK127","CK129","CK135","CK137","CK141"],
    'hf_ckd' : ["CK116","CK126","CK136","CK138"]
}

# Initialize sumary variables
total_df = pd.DataFrame(columns=['n_genes_by_counts','total_counts',
                                   'pct_counts_mt','doublet_score'])
total_subsets = np.zeros(7).astype(np.int)
summary_df = []

# Run QC plots for each sample and store summary
for condition, samples in meta.items():
    for sample in samples:
        print(sample)
        input_path = '../plot_data/qc/{0}.pkl'.format(sample)
        plot_data = pickle.load(open(input_path, "rb"))
        
        # Filter params
        mt_thr = plot_data['mt_thr']
        gene_qnt = plot_data['gene_qnt']
        doublet_thr = plot_data['doublet_thr']
        df = plot_data['df']
        gene_thr = np.quantile(df.n_genes_by_counts, gene_qnt)
        
        # Create QC plots
        fig, axes = plt.subplots(2,2, figsize=(8,6))
        fig.suptitle('QC metrics {0}'.format(sample), fontsize=11)
        axes = axes.flatten()
        
        # Plot MT
        plot_mt_vs_counts(df, axes[0], mt_thr=mt_thr)
        
        # Plot ngenes
        plot_ngenes_vs_counts(df, axes[1], gene_thr=gene_thr)
        
        # Plot doublet score
        plot_doublet_scores(df, axes[2], doublet_thr=doublet_thr)
        
        # Plot filter overlap
        subsets = compute_overlaps(df, mt_thr=mt_thr, gene_thr=gene_thr, doublet_thr=doublet_thr)
        plot_venn(subsets, axes[3])
        
        # Adjust plots
        fig.tight_layout()
        fig.subplots_adjust(top=0.88)
        fig.set_facecolor('white')
        
        # Write to png
        fig.savefig('../plots/qc_{0}.png'.format(sample))
        
        # Filter and append
        msk = (df.n_genes_by_counts < gene_thr) & \
          (df.pct_counts_mt < mt_thr) &  \
          (df.doublet_score < doublet_thr)
        df = df[msk]
        
        total_df = total_df.append(df, ignore_index=True)
        summary_df.append([sample, len(msk), sum(msk)])
        
        # Update total subsets
        total_subsets = total_subsets + subsets


# Summary plots
fig = plt.figure(figsize=(8, 6))
fig.suptitle('Merged QC metrics', fontsize=11)
gs = fig.add_gridspec(2, 3)

ax = fig.add_subplot(gs[0,0])
plot_mt_vs_counts(total_df, ax, mt_thr=np.nan)

ax = fig.add_subplot(gs[0,1])
plot_ngenes_vs_counts(total_df, ax, gene_thr=np.nan)

ax = fig.add_subplot(gs[0,2])
plot_doublet_scores(total_df, ax, doublet_thr=np.nan)

ax = fig.add_subplot(gs[1,0:2])
pd.DataFrame(summary_df, columns=["Sample ID", "Before", "After"]).plot.bar(x=0, ax=ax)

ax = fig.add_subplot(gs[1,2])
plot_venn(total_subsets, ax)

fig.tight_layout()
fig.subplots_adjust(top=0.88)
fig.set_facecolor('white')

# Save
fig.savefig('../plots/qc_summary.png')
