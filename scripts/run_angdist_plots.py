import pandas as pd
import numpy as np

import pickle

import matplotlib.pyplot as plt
import seaborn as sns

plot_data = pickle.load(open("../plot_data/angdist/plot_data.pkl", "rb"))

fig, axes = plt.subplots(2,3, figsize=(12,9), sharey=True, dpi=150)
fig.suptitle('Transcriptomical shifts', fontsize=15)
axes = axes.flatten()
i = 0
for k,df in plot_data.items():
    ax = axes[i]
    order = df.groupby('cell_type').median().sort_values('dists', ascending=False).index
    sns.boxplot(data=df, x='cell_type', y='dists', ax=ax, order=order, palette='coolwarm_r',linewidth=0.5)
    ax.axhline(1, ls='--', c='black')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    ax.set_title(k, fontsize=11)
    ax.set_xlabel('')
    ax.set_ylabel('')
    i +=1

fig.subplots_adjust(right=0.9)
fig.tight_layout()
fig.set_facecolor('white')

# Save
fig.savefig('../plots/angdist.png')