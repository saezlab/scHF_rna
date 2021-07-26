import pandas as pd
import numpy as np

import scanpy as sc

from utils import get_design, get_contrast, limma_fit


# Read data
input_path = '../qc_data/bulk.h5ad'
adata = sc.read_h5ad(input_path).raw.to_adata()
# Need to set to str othwerwise pivot crashes
adata.obs['condition'] = adata.obs['condition'].astype(str)
data = pd.DataFrame(adata.X.T, index=adata.var.index, columns=adata.obs.index)
data.columns.name = None

### Contrasts to test ###

contr_dict = {
    'hf-healthy' : ['hf','healthy'],
    'hf_ckd-healthy' : ['hf_ckd','healthy'],
    'hf_ckd-hf' : ['hf_ckd','hf']
}

###

# DEG
design = get_design(adata.obs, 'condition')
contr_matrix = get_contrast(design, contr_dict)
data = pd.DataFrame(adata.X.T, index=adata.var.index, columns=adata.obs.index)
design = design.loc[data.columns]
deg = limma_fit(data, design, contr_matrix).sort_values(['contrast', 'pvals'])

# Save
deg.to_csv('../plot_data/deg/bulk_deg.csv', index=False)