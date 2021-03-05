import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

# Read DEG df
input_path = '../plot_data/deg/deg.csv'
df = pd.read_csv(input_path)

def plot_volcano(name, df, ax, max_num=np.max(np.abs(df['logfoldchanges'])), s=10, fontsize=12):
    # Get lfc and pvals
    lfc = df['logfoldchanges']
    pvals = -np.log10(df['pvals'])
    
    # Mask significant genes
    msk = (pvals > -np.log10(0.05)) & (np.abs(lfc) > 0.5)
    
    # Plot scatter
    ax.set_title(name)
    ax.scatter(lfc[~msk], pvals[~msk], c='gray', s=s)
    ax.scatter(lfc[msk], pvals[msk], c='red', s=s)
    ax.set_xlim(-max_num, max_num)
    ax.set_xlabel('Logfoldchanges', fontsize=fontsize)
    ax.set_ylabel('-log10(pvalue)', fontsize=fontsize)
    ax.set_box_aspect(1)


conditions = np.unique(df['condition'])
cell_types = np.unique(df['cell_type'])

for condition in conditions:
    # Define figure
    fig, axes = plt.subplots(3,4, figsize=(4*3, 3*3) , dpi=150, sharey=True, sharex=True)
    axes = axes.flatten()
    fig.suptitle('{0} vs healthy'.format(condition), fontsize=16)
    axes[0].set_visible(False) 
    
    # Get max lfc to set as limit for plot
    max_lfc = np.max(np.abs(df[df['condition']==condition]['logfoldchanges']))
    
    # Volcano for each cell type
    for ctype,ax in zip(cell_types, axes[1:]):
        deg = df[(df['condition'] == condition) & (df['cell_type'] == ctype)]
        if deg.shape[0] == 0:
            ax.set_visible(False) 
            continue
        plot_volcano(ctype, deg, ax, max_lfc)
    
    # Save figure
    fig.tight_layout()
    fig.set_facecolor('white')
    fig.savefig('../plots/deg_{0}'.format(condition))

axes[0].set_visible(False) 
    
fig.tight_layout()
fig.set_facecolor('white')