import scanpy as sc
import numpy as np
import pandas as pd

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat

from utils import get_count_df


# Define order of categories and reference cell type
order_cat = ['healthy', 'hf', 'hf_ckd', 'acidosis']
reference_cell_type = 'neuronal'

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Get cell type count df
df = get_count_df(adata, drop=['unknown'])

# Create scCODA AnnData object
data = dat.from_pandas(df, covariate_columns=['sample_id', 'condition'])
data.obs['condition'] = data.obs['condition'].astype('category').cat.reorder_categories(order_cat)

# Initialize empty df
eff_df = pd.DataFrame(columns=data.var.index)

# Test effects for each condition that is not healthy
for cond in order_cat[1:]:
    # Select control and condition data
    data_cond = data[data.obs["condition"].isin(["healthy", cond])]
    model_cond = mod.CompositionalAnalysis(data_cond, formula="condition", reference_cell_type=reference_cell_type)
    # Run MCMC
    sim_results = model_cond.sample_hmc()
    # Add to df
    eff_df.loc[cond] = sim_results.summary_prepare()[1].reset_index(level=[0])['Final Parameter']

# Add effects to adata object
adata.uns['coda'] = eff_df.T


# Write to file
adata.write('../qc_data/integrated.h5ad')