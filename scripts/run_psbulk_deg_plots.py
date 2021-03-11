import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

# Read DEG df
input_path = '../plot_data/deg/deg.csv'
df = pd.read_csv(input_path)

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
    fig.savefig('../plots/psbulk_deg_{0}'.format(condition))

axes[0].set_visible(False) 
    
fig.tight_layout()
fig.set_facecolor('white')