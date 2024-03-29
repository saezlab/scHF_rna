import scanpy as sc
import numpy as np
import pandas as pd

from sccoda.util import cell_composition_data as dat
from plotting import stacked_barplot

import matplotlib.pyplot as plt
import seaborn as sns

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

# Get scCODA input:
df = adata.obs.pivot(index=['sample_id','condition'], columns='cell_type', values='cell_num')
df[np.isnan(df)] = 0
df = df.rename(columns=str).reset_index()
df.columns.name = None
data = dat.from_pandas(df, covariate_columns=['sample_id', 'condition'])

# Load eff df
eff_df = adata.uns['coda']

# Plot proportions and CODA effects
fig, ax = plt.subplots(1,2,figsize=(10,5), dpi=100)
ax=ax.flatten()
stacked_barplot(data, feature_name="condition", ax=ax[0])
sns.heatmap(eff_df, center=0, cmap='coolwarm', annot=True, square=True, ax=ax[1])

# Adjust
fig.tight_layout()
fig.set_facecolor('white')

# Save
fig.savefig('../plots/coda.png')