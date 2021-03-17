import scanpy as sc
import numpy as np
import pandas as pd

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat


# Define reference cell type
reference_cell_type = 'neuronal'

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

# Get unique conditions
conditions = np.unique(adata.obs['condition'])

# Get scCODA input:
# sample_id  condition  T-cells  adipocytes  cardiomyocyte
# CK114      healthy    238      0           1246.0
df = adata.obs.pivot(index=['sample_id','condition'], columns='cell_type', values='cell_num')
df[np.isnan(df)] = 0
df = df.rename(columns=str).reset_index()
df.columns.name = None

# Create scCODA AnnData object
data = dat.from_pandas(df, covariate_columns=['sample_id', 'condition'])

# Initialize empty df
eff_df = pd.DataFrame(columns=data.var.index)

# Test effects for each condition that is not healthy
for i,cond_a in enumerate(conditions):
    for j,cond_b in enumerate(conditions):
        if i < j:
            contrast = '{0}-{1}'.format(cond_b, cond_a)
            print(contrast)
            
            # Select control and condition data
            data_cond = data[data.obs["condition"].isin([cond_a, cond_b])]
            data_cond.obs['condition'] = data_cond.obs['condition'].astype(str)
            data_cond.obs['condition'] = data_cond.obs['condition'].astype('category').cat.reorder_categories([cond_a, cond_b])
            model_cond = mod.CompositionalAnalysis(data_cond, formula="condition", reference_cell_type=reference_cell_type)
            
            # Run MCMC
            sim_results = model_cond.sample_hmc()
            
            # Add to df
            eff_df.loc[contrast] = sim_results.summary_prepare()[1].reset_index(level=[0])['Final Parameter']

# Add effects to adata object
adata.uns['coda'] = eff_df.T


# Write to file
adata.write('../qc_data/pseudobulk.h5ad')