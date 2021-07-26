import pandas as pd
import numpy as np

import pickle
from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns

plot_data = pickle.load(open("../plot_data/angdist/plot_data.pkl", "rb"))

fig, axes = plt.subplots(1,3, figsize=(12,4), sharey=True, dpi=150)
fig.suptitle('Transcriptomic shifts', fontsize=15)
axes = axes.flatten()
i = 0

for k,df in plot_data.items():
    ax = axes[i]
    order = df.groupby('cell_type').median().sort_values('dists', ascending=False).index
    
    # Boxplots
    sns.boxplot(data=df, x='cell_type', y='dists', ax=ax, order=order, palette='coolwarm_r',linewidth=0.5)
    ax.axhline(1, ls='--', c='black')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_title(k, fontsize=11)
    ax.set_xlabel('')
    ax.set_ylabel('inter dist / intra dist')
    
    # P-vals
    for j, cell_type in enumerate(order):
        dists = df[df['cell_type'] == cell_type]['dists'].values
        _, pval = stats.ttest_1samp(dists, 1, alternative='greater')
        if pval < 0.05:
            pval = '*'
        else:
            pval = ''
        
        ax.text(j, 3, pval)
    
    i +=1

fig.subplots_adjust(right=0.9)
fig.tight_layout()
fig.set_facecolor('white')

# Save
fig.savefig('../plots/psbulk_angdist.png')