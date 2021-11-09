import scanpy as sc
import numpy as np
import pandas as pd

from anndata import AnnData
from utils import vsn_normalize

'''Generates a new AnnData object of pseudobulk expression profiles per cell type and sample.'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
meta = sc.read_h5ad(input_path)

# Get complete list of genes
g_arr = meta.raw.var.index.values

# Filter out unknown cells
meta = meta[meta.obs['cell_type'] != 'unknown']

# Get metadata
meta = meta.obs

# Get cell types
cell_types = np.unique(meta['cell_type'])

# Get samples ids
sample_ids = np.unique(meta['sample_id'])

# Get min size of groups
min_size = meta.drop_duplicates(["sample_id", "condition"])['condition'].value_counts().min()

lst_msks = []
X = []
col_sample_id = []
col_condition = []
col_cell_type = []
col_cell_num = []
col_cell_prop = []

for sample_id in sample_ids:
    # Get condition
    condition = meta[meta['sample_id']==sample_id]['condition'].head(1).values[0]
    
    # Open RNA counts and filter by overlaped cells and genes
    input_path = '../data/{0}/outs/filtered_feature_bc_matrix/'.format(sample_id)
    raw_adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', cache=True)
    raw_adata = raw_adata[:, g_arr]

    # Subset meta
    s_meta = meta[meta['sample_id']==sample_id]
    
    for cell_type in cell_types:
        # Get cell ids
        cell_ids = s_meta[s_meta['cell_type']==cell_type].index
        cell_ids = ['-'.join(c_id.split('-')[:-1]) for c_id in cell_ids]

        # If sample has less than 10 cells skip
        if len(cell_ids) >= 10:
            # Subset
            cell_adata = raw_adata[cell_ids]
        else:
            continue
        
        # Get 5% msk
        expr_msk = cell_adata.X.A > 0
        expr_msk = np.sum(expr_msk, axis=0) / expr_msk.shape[0]
        expr_msk = expr_msk > 0.05
        
        # Get psbulk profile
        expr_X = np.sum(cell_adata.X.A, axis=0)
        
        # If sample ahs less than 1000 reads skip
        if np.sum(expr_X) < 1000:
            continue

        # Append
        lst_msks.append(expr_msk)
        X.append(expr_X)
        col_sample_id.append(sample_id)
        col_condition.append(condition)
        col_cell_type.append(cell_type)
        col_cell_num.append(cell_adata.shape[0])
        col_cell_prop.append(cell_adata.shape[0]/raw_adata.shape[0])

lst_msks = np.array(lst_msks)

# Crete AnnData object
pb_adata = AnnData(np.array(X))
pb_adata.obs.index = [sample_id+'_'+cell_type for sample_id,cell_type in zip(col_sample_id, col_cell_type)]
pb_adata.obs['sample_id'] = col_sample_id
pb_adata.obs['condition'] = col_condition
pb_adata.obs['cell_type'] = col_cell_type
pb_adata.obs['cell_num'] = col_cell_num
pb_adata.obs['cell_prop'] = col_cell_prop
pb_adata.var.index = g_arr

# Remove genes that are less than 5% in min_size samples
for cell_type in cell_types:
    # Get cell type and gene msks
    ctype_msk = pb_adata.obs['cell_type'] == cell_type
    msk = np.sum(lst_msks[ctype_msk], axis=0) >= min_size
    print(cell_type, 'num genes:', np.sum(msk))
    
    # Subset
    sub_X = np.array(pb_adata[ctype_msk, msk].X)
    
    # VSN normalize
    X_norm = vsn_normalize(sub_X)
    
    # Set filtered genes expr to 0
    Xdim = (np.sum(ctype_msk), msk.shape[0])
    sub_X = np.zeros(Xdim)
    
    # Add normalized values to matrix
    sub_X[:,msk] = X_norm
    
    # Update expr matrix in Anndata object
    pb_adata[ctype_msk].X = sub_X
    
# Write to file
pb_adata.write('../qc_data/pseudobulk.h5ad')
