import numpy as np
import pandas as pd

import pickle

import matplotlib.pyplot as plt
import seaborn as sns


input_path = '../plot_data/jsd/plot_data.pkl'
plot_data = pickle.load(open(input_path,"rb"))
samples_ids = plot_data['samples_ids']

conditions = plot_data['conditions']
cell_types = plot_data['cell_types']
jsds = plot_data['jsds']
coords = plot_data['coords']

fig, axes = plt.subplots(3, 4, figsize=(12,9), dpi=150)
fig.suptitle('JSD distances MDS', fontsize=16)
axes = axes.flatten()

for i in range(len(cell_types)):
    ax = axes[i]
    ctype = cell_types[i]
    ax.set_title(ctype)
    max_num = np.nanmax(np.abs(coords[i]))
    max_num = max_num + max_num * 0.1
    for cond in np.unique(conditions):
        msk = conditions == cond
        coord = coords[i][msk]
        ax.scatter(coord[:,0], coord[:,1], label=cond, s=20)
        ax.set_box_aspect(1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-max_num,max_num)
        ax.set_ylim(-max_num,max_num)
        ax.set_xlabel('MDS 1')
        ax.set_ylabel('MDS 2')
        
fig.tight_layout()
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='right', frameon=False, fontsize=11)
fig.subplots_adjust(right=0.9)
fig.set_facecolor('white')

fig.savefig('../plots/psbulk_jsd_mds.png')

fig, axes = plt.subplots(3, 4, figsize=(12,9), dpi=150, sharex=True, sharey=True)
fig.suptitle('JSD distances MDS', fontsize=16)
axes = axes.flatten()

for i in range(len(cell_types)):
    ax = axes[i]
    ctype = cell_types[i]
    ax.set_title(ctype)
    jsd = jsds[i]
    order = np.argsort(conditions)
    sns.heatmap(jsd[order][:,order], ax=ax, square=True, xticklabels=samples_ids[order], yticklabels=samples_ids[order],
               cbar_kws={"shrink": .5}, cmap='viridis', robust=True)
        
fig.tight_layout()
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='right', frameon=False, fontsize=11)
fig.set_facecolor('white')

fig.savefig('../plots/psbulk_jsd_dist.png')