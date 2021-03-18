import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm

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

    
def stacked_barplot(data, feature_name, ax, cmap=cm.tab20):
    # cell type names
    type_names = data.var.index
    levels = pd.unique(data.obs[feature_name])
    n_levels = len(levels)
    feature_totals = np.zeros([n_levels, data.X.shape[1]])

    for level in range(n_levels):
        l_indices = np.where(data.obs[feature_name] == levels[level])
        feature_totals[level] = np.sum(data.X[l_indices], axis=0)
        
    y = feature_totals
    title=feature_name
    level_names=levels
    
    n_bars, n_types = y.shape
    r = np.array(range(n_bars))
    sample_sums = np.sum(y, axis=1)

    barwidth = 0.85
    cum_bars = np.zeros(n_bars)

    for n in range(n_types):
        bars = [i / j * 100 for i, j in zip([y[k][n] for k in range(n_bars)], sample_sums)]
        ax.bar(r, bars, bottom=cum_bars, color=cmap(n % cmap.N), width=barwidth, label=type_names[n], linewidth=0)
        cum_bars += bars
    ax.set_title(title)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    ax.set_xticks(r)
    ax.set_xticklabels(level_names, rotation=45)
    ax.set_ylabel("Proportion")

def volcano(name, lfc, pvals, ax, max_num=None,
                 lfc_thr=0.5, p_thr=0.05, s=10, fontsize=12):
    '''Volcano plot from a list of lfc and untrasnformed pvalues'''

    if max_num is None:
        max_num=np.max(np.abs(lfc))
    
    # Transform pvals
    pvals = -np.log10(pvals)
    
    # Mask significant genes
    msk = (pvals > -np.log10(p_thr)) & (np.abs(lfc) > lfc_thr)
    
    # Plot scatter
    ax.set_title(name)
    ax.scatter(lfc[~msk], pvals[~msk], c='gray', s=s)
    ax.scatter(lfc[msk], pvals[msk], c='red', s=s)
    ax.set_xlim(-max_num, max_num)
    ax.set_xlabel('LogFC', fontsize=fontsize)
    ax.set_ylabel('-log10(pvalue)', fontsize=fontsize)
    ax.set_box_aspect(1)
    
def dotplot(ref, df, num=30, fontsize=9, figsize=(12,6)):
    # Define figure
    fig, ax = plt.subplots(1,1, dpi=150, figsize=figsize)
    ax.set_title(ref, fontsize=fontsize+5)
    
    # Add grid and set it to background
    ax.grid(True)
    ax.set_axisbelow(True)
    
    # Dot plot
    max_num = np.max(np.abs(df['coeff']))
    sc = ax.scatter(
        x=df['name'],
        y=df['ctype.cond'],
        c=df['coeff'],
        s=-np.log(df['pvals']) * num,
        cmap='coolwarm',
        vmax=max_num,
        vmin=-max_num,
        edgecolor=['black' if pval < 0.05 else 'white' for pval in df['pvals']]
    )
    
    # Format dot plot ticks
    ax.tick_params(axis='x', rotation=90, labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)

    # Plot pvalue dot sizes legend
    handles, labels = sc.legend_elements("sizes", num=4)
    labels = ['$\\mathdefault{'+'{0:.2f}'.format(int(label.split('{')[1].split('}')[0])/num)+'}$' for label in labels]
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1), frameon=False, title='-log(pvalue)')
    
    # Add color bar
    cax = fig.add_axes([0.945, 0.25, 0.025, 0.35])
    cbar = fig.colorbar(sc, cax=cax, orientation='vertical')
    cbar.ax.set_title('Mean change')
    
    # Format figure
    fig.tight_layout()
    fig.set_facecolor('white')
    
    return fig

def plot_ora(name, df, ax, top=10, fontsize=11):
    df = df.sort_values('adj_pvalues', ascending=True).head(top)
    names = np.flip(df['descr'].tolist())
    pvals = np.flip(-np.log10(df['adj_pvalues']))
    ax.barh(names, pvals, color='gray')
    ax.axvline(x=-np.log10(0.05), c='black', ls='--')
    ax.set_xlabel('-log10(adj_pval)', fontsize=fontsize)
    ax.set_title(name, fontsize=fontsize)