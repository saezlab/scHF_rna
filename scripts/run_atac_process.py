import scanpy as sc
import numpy as np
import pandas as pd

from anndata import AnnData

# Read meta and ATAC peaks
meta = pd.read_csv('../data/metadata.csv').set_index('code')
adata = AnnData(pd.read_csv('../data/ATAC/psbulk_peaks.csv', index_col=0).T)

# Filter and transform
sc.pp.filter_genes(adata, min_cells=3)
#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)

# Transform metadata to same format
conditions = []
cell_types = []
sample_ids = []
for sample in adata.obs.index:
    sample = sample.split('_')
    if len(sample) == 2:
        code, cell_type = sample
    else:
        code = sample[0]
        cell_type = '_'.join(sample[1:])
    cell_types.append(cell_type)
    condition = meta.loc[code, 'condition']
    conditions.append(condition)
    sample_id = meta.loc[code, 'sample_id']
    sample_ids.append(sample_id)
adata.obs['condition'] = conditions
adata.obs['cell_type'] = cell_types
adata.obs['sample_id'] = sample_ids


# Rewrite indexes to follow same format
adata.obs['index'] = adata.obs['sample_id'] + '_' + adata.obs['cell_type']
adata.obs = adata.obs.reset_index().set_index('index')

# Write
adata.write('../qc_data/atac_peaks.h5ad')
