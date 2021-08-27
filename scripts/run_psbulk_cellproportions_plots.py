import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from plotting import cond_colors, hextofloats

df = pd.read_csv('../plot_data/var_shifts/cell_proportions.csv')
table = pd.read_csv('../tables/cellproportions.csv')
table = table.set_index(['condition','cell_type'])

order = np.sort(df['CellType'].unique())
hue_order = list(cond_colors.keys())



fig, ax = plt.subplots(1,1, figsize=(9,3), tight_layout=True, dpi=150, facecolor='white')
sns.boxplot(x="CellType", y="Prop", hue="Condition", fliersize=0, palette=cond_colors,
            data=df, ax=ax, linewidth=0.5, order=sorted(df['CellType'].unique()))
ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
ax.set_title('Cell type composition shifts', fontsize=11)
ax.set_xlabel('')
ax.set_ylabel('Proportion')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

n = len(cond_colors)
cell_types = np.repeat(order, n)
for i,patch in enumerate(ax.artists):
    hue = hue_order[i % n]
    ctype = cell_types[i]
    adj_pval = table.loc[hue].loc[ctype]['adj_pval']
    if adj_pval < 0.05:
        a = 1
    else:
        a = 0.33
    r, g, b = hextofloats(cond_colors[hue])
    patch.set_facecolor((r, g, b, a))

fig.savefig('../plots/psbulk_cellproportions.pdf', bbox_inches='tight')