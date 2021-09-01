import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from plotting import cond_colors 

# Load data
df = pd.read_csv('../plot_data/jsd/jsd.csv')

types = np.unique(df['type'])
cell_types = np.unique(df['cell_type'])

for typ in types:
    data = df[df['type'] == typ]

    fig, axes = plt.subplots(3,4, figsize=(12,9), tight_layout=True, facecolor='white')
    axes = axes.flatten()

    for i,cell_type in enumerate(cell_types):
        ax = axes[i]
        subadata = data[data.cell_type == cell_type]
        sns.scatterplot(data=subadata, ax=ax,
                        x='x', y='y', hue='condition', palette=cond_colors, legend=False)
        ax.set_box_aspect(1)
        ax.set_title(cell_type)
        if not np.all(np.isnan(subadata[['x','y']].values)):
            n_max = np.nanmax(np.abs(subadata[['x','y']].values))
            n_max = n_max + n_max * 0.1
            ax.set_xlim(-n_max,n_max)
            ax.set_ylim(-n_max,n_max)
        ax.axvline(x=0, alpha=0.25, color='grey', linestyle='--', zorder=0)
        ax.axhline(y=0, alpha=0.25, color='grey', linestyle='--', zorder=0)
    fig.savefig('../plots/psbulk_jsd_dist_{0}.pdf'.format(typ))