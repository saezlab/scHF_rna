import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

from plotting import cond_colors 

# Load data
jsds = pd.read_csv('../plot_data/jsd/jsd.csv')

types = np.unique(jsds['type'])
cell_types = np.unique(jsds['cell_type'])

dsts = pd.read_csv('../plot_data/var_shifts/shifts.csv')

dpvl = pd.read_csv('../tables/shifts_celltypes.csv')
dpvl = dpvl.set_index(['condition','cell_type'])

prps = pd.read_csv('../plot_data/var_shifts/cell_proportions.csv')
ppvl = pd.read_csv('../tables/cellproportions.csv')
ppvl = ppvl.set_index(['condition','cell_type'])

fig = plt.figure(constrained_layout=True, figsize=(12,5), facecolor='white')
gs = fig.add_gridspec(2, 5)

order = list(cond_colors.keys())


ax = fig.add_subplot(gs[0:2,0:2])
cell_type = 'All'
data = jsds[(jsds['type'] == 'BOTH')&(jsds['cell_type'] == cell_type)]
sns.scatterplot(data=data, ax=ax, x='x', y='y', hue='condition', palette=cond_colors, legend=True, s=1)
for i in range(data.shape[0]):
    cond = data.condition.values[i]
    ax.text(x=data.x.values[i],
            y=data.y.values[i],
            s=data.sample_id.values[i],
            c=cond_colors[cond]
           )
ax.set_title(cell_type)
n_max = np.nanmax(np.abs(data[['x','y']].values))
n_max = n_max + n_max * 0.1
ax.set_xlim(-n_max,n_max)
ax.set_ylim(-n_max,n_max)
ax.axvline(x=0, alpha=0.25, color='grey', linestyle='--', zorder=0)
ax.axhline(y=0, alpha=0.25, color='grey', linestyle='--', zorder=0)
ax.set_title('')
ax.set_xlabel('MDS1')
ax.set_ylabel('MDS2')
ax.legend(frameon=False)

# Cardios
cell_type = 'cardiomyocyte'
ax = fig.add_subplot(gs[0,2])
data = prps[prps['CellType'] == cell_type]
sns.boxplot(x="Condition", y="Prop", fliersize=0, palette=cond_colors,
            data=data, ax=ax, linewidth=1, order=order)
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.set_title(cell_type, fontsize=11)
ax.set_xlabel('')
ax.set_ylabel('Proportions')
for i,cond in enumerate(cond_colors.keys()):
    adj_pval = ppvl.loc[cond].loc[cell_type]['adj_pval']
    if adj_pval < 0.05:
        ax.artists[i].set_hatch('\\\\')
ax = fig.add_subplot(gs[1,2])
ax.set_box_aspect(1)
data = dsts[dsts['cell_type'] == cell_type]
sns.boxplot(data=data, x='condition', y='dist', ax=ax, palette=cond_colors, linewidth=1, 
                fliersize=0, order=order)
ax.axhline(1, ls='--', c='black')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.set_xlabel('')
ax.set_ylabel('Molecular distances')
for i,cond in enumerate(cond_colors.keys()):
    adj_pval = dpvl.loc[cond].loc[cell_type]['adj_pval']
    if adj_pval < 0.05:
        ax.artists[i].set_hatch('\\\\')

# Fibros
cell_type = 'fibroblast'
ax = fig.add_subplot(gs[0,3])
data = prps[prps['CellType'] == cell_type]
sns.boxplot(x="Condition", y="Prop", fliersize=0, palette=cond_colors,
            data=data, ax=ax, linewidth=1, order=order)
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.set_title(cell_type, fontsize=11)
ax.set_xlabel('')
ax.set_ylabel('')
for i,cond in enumerate(cond_colors.keys()):
    adj_pval = ppvl.loc[cond].loc[cell_type]['adj_pval']
    if adj_pval < 0.05:
        ax.artists[i].set_hatch('\\\\')
ax = fig.add_subplot(gs[1,3])
ax.set_box_aspect(1)
data = dsts[dsts['cell_type'] == cell_type]
sns.boxplot(data=data, x='condition', y='dist', ax=ax, palette=cond_colors, linewidth=1, 
                fliersize=0, order=order)
ax.axhline(1, ls='--', c='black')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.set_xlabel('')
ax.set_ylabel('')
for i,cond in enumerate(cond_colors.keys()):
    adj_pval = dpvl.loc[cond].loc[cell_type]['adj_pval']
    if adj_pval < 0.05:
        ax.artists[i].set_hatch('\\\\') 
         
# Endos
cell_type = 'endothelial'
ax = fig.add_subplot(gs[0,4])
data = prps[prps['CellType'] == cell_type]
sns.boxplot(x="Condition", y="Prop", fliersize=0, palette=cond_colors,
            data=data, ax=ax, linewidth=1, order=order)
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.set_title(cell_type, fontsize=11)
ax.set_xlabel('')
ax.set_ylabel('')
for i,cond in enumerate(cond_colors.keys()):
    adj_pval = ppvl.loc[cond].loc[cell_type]['adj_pval']
    if adj_pval < 0.05:
        ax.artists[i].set_hatch('\\\\')
ax = fig.add_subplot(gs[1,4])
ax.set_box_aspect(1)
data = dsts[dsts['cell_type'] == cell_type]
sns.boxplot(data=data, x='condition', y='dist', ax=ax, palette=cond_colors, linewidth=1, 
                fliersize=0, order=order)
ax.axhline(1, ls='--', c='black')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.set_xlabel('')
ax.set_ylabel('')
for i,cond in enumerate(cond_colors.keys()):
    adj_pval = dpvl.loc[cond].loc[cell_type]['adj_pval']
    if adj_pval < 0.05:
        ax.artists[i].set_hatch('\\\\')
        
fig.savefig('../plots/mds_prop_shifts.pdf')