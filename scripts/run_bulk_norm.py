import pandas as pd
import numpy as np

import scanpy as sc
from anndata import AnnData

from utils import vsn_normalize

### Samples to filter ###

samples_filter = [] #['CAD10', 'N15', 'RF17']

###

# Read data
meta = pd.read_csv('../data/metadata.csv')
meta = meta[meta.bulk_id != 'None']
meta.index = meta.bulk_id.values
adata = pd.read_csv('../data/GSE160145.csv', index_col=0)
adata = adata[meta.index]

# Create AnnData object
adata = AnnData(adata.values.T, var=pd.DataFrame(index=adata.index), obs=meta)

# General filter
sc.pp.filter_cells(adata, min_genes=5000)
sc.pp.filter_genes(adata, min_cells=3)

# Filter samples
if len(samples_filter) > 0:
    msk = np.array([s_id not in samples_filter for s_id in adata.obs.index])
    adata = adata[msk]

# Filter lowly expressed genes
adata.X[adata.X == 0] = np.nan
adata.X[np.log2(adata.X) < 1] = np.nan
msk = np.sum(np.isnan(adata.X), axis=0) < 2
adata = adata[:,msk]
adata.X[np.isnan(adata.X)] = 0
adata.X = adata.X.astype(int)

# VSN normalize and store
adata.X = vsn_normalize(adata.X)
adata.raw = adata

# Compute PCA and UMAP
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=adata.shape[0]-1)
sc.tl.umap(adata)

# Save
adata.write('../qc_data/bulk.h5ad')