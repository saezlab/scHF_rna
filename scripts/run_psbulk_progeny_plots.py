import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


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
    labels = ['$\\mathdefault{'+str(int(label.split('{')[1].split('}')[0])/num)+'}$' for label in labels]
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1), frameon=False, title='-log(pvalue)')
    
    # Add color bar
    cax = fig.add_axes([0.945, 0.25, 0.025, 0.35])
    cbar = fig.colorbar(sc, cax=cax, orientation='vertical')
    cbar.ax.set_title('Mean change')
    
    # Format figure
    fig.tight_layout()
    fig.set_facecolor('white')
    
    return fig


# Read teste pathways
results = pd.read_csv('../plot_data/func/progeny.csv')
results['ctype.cond'] = results['cell_type'] + ' ' + results['cond']
refs = np.unique(results['ref'])

for ref in refs:
    # Filter by a reference condition
    df = results[(results['ref'] == ref)]
    
    # Get significant elements
    sign = df[df['pvals'] < 0.05]
    names = np.unique(sign['name'].tolist())
    ctconds = np.unique(sign['ctype.cond'].tolist())
    
    # Filter by combination of significant cell type and condition
    msk_ctconds = np.array([ctcond in ctconds for ctcond in df['ctype.cond'].tolist()])
    msk_names = np.array([name in names for name in df['name'].tolist()])
    df = df[msk_ctconds * msk_names]
    msk = df['pvals'] > 0.05
    
    # Hide non significant
    df['pvals'].loc[msk] = 1
    df['coeff'].loc[msk] = 0
    
    # Plot dot plot
    fig = dotplot(ref, df)
    fig.savefig('../plots/psbulk_progeny_{0}.png'.format(ref), bbox_inches='tight')