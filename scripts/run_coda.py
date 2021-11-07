import scanpy as sc
import numpy as np
import pandas as pd

import os

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat


# Define reference cell type
reference_cell_type = 'neuronal'

# Get scCODA input:
# sample_id  condition  T-cells  adipocytes  cardiomyocyte
# CK114      healthy    238      0           1246.0
df = pd.read_csv('../plot_data/var_shifts/cell_proportions.csv')
df = df.pivot(index=['Patient', 'Data', 'Condition'], columns='CellType', values='Counts')
df[np.isnan(df)] = 0.5
df = df.rename(columns=str).reset_index()
df.columns.name = None
del df['Data']

# Create scCODA AnnData object
data = dat.from_pandas(df, covariate_columns=['Patient', 'Condition'])

# Get unique contrasts
conditions = np.unique(df['Condition'])

# Initialize empty df
eff_df = pd.DataFrame(columns=data.var.index)

# Contrasts
contrasts = ['HF-AvsHealthy', 'HF-CKDvsHealthy','HF-CKDvsHF-A']

# Test effects for each condition that is not healthy
for contrast in contrasts:
    print(contrast)
    cond_a, cond_b = contrast.split('vs')

    # Select control and condition data
    data_cond = data[data.obs["Condition"].isin([cond_b, cond_a])]
    data_cond.obs['Condition'] = data_cond.obs['Condition'].astype(str)
    data_cond.obs['Condition'] = data_cond.obs['Condition'].astype('category').cat.reorder_categories([cond_b, cond_a])
    model_cond = mod.CompositionalAnalysis(data_cond, formula="Condition", reference_cell_type=reference_cell_type)

    # Run MCMC
    sim_results = model_cond.sample_hmc()

    # Add to df
    eff_df.loc[contrast] = sim_results.summary_prepare()[1].reset_index(level=[0])['Final Parameter']

# Save
os.makedirs('../plot_data/coda', exist_ok=True)
eff_df.to_csv('../plot_data/coda/coda.csv')
