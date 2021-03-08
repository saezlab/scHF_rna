import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


def plot_ora(name, df, ax, top=10, fontsize=11):
    df = df.sort_values('adj_pvalues', ascending=True).head(top)
    names = np.flip(df['descr'].tolist())
    pvals = np.flip(-np.log10(df['adj_pvalues']))
    ax.barh(names, pvals, color='gray')
    ax.axvline(x=-np.log10(0.05), c='black', ls='--')
    ax.set_xlabel('-log10(adj_pval)', fontsize=fontsize)
    ax.set_title(name, fontsize=fontsize)
    
    
# Read ORA df
input_path = '../plot_data/deg/ora.csv'
df = pd.read_csv(input_path)

# Get different conditions and cell types
conditions = np.unique(df['condition'])

for condition in conditions:
    ora = df[df['condition'] == condition]
    cell_types = np.unique(ora['cell_type'])
    fig, axes = plt.subplots(len(cell_types),1, dpi=150, figsize=(4, 4*len(cell_types)), sharex=True)
    fig.suptitle('ORA DEG {0}.vs.healty'.format(condition), fontsize=11)
    axes = axes.flatten()
    for cell_type, ax in zip(cell_types, axes):
        plot_ora(cell_type, ora[ora['cell_type'] == cell_type], ax=ax)
    
    fig.subplots_adjust(top=0.95)
    fig.tight_layout()
    fig.set_facecolor('white')
    
    fig.savefig('../plots/ora_{0}.png'.format(condition), bbox_inches='tight')