import scanpy as sc
import numpy as np
import pandas as pd

from anndata import AnnData

'''Generates a new AnnData object of pseudobulk expression profiles per cell type and sample.'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)
# Filter out unknown cells
adata = adata[adata.obs['cell_type'] != 'unknown']


#adata.obs['condition'] = adata.obs['condition'].astype('str')
#adata.obs['condition'][adata.obs['condition'] != 'hf_ckd'] = "other"
#adata.obs['condition'] = adata.obs['condition'].astype('category').cat.reorder_categories(['other', 'hf_ckd'])

X = []
col_sample_id = []
col_condition = []
col_cell_type = []
col_cell_num = []
col_cell_prop = []

# Get samples ids
sample_ids = np.unique(adata.obs['sample_id'])
for sample_id in sample_ids:
    print(sample_id)
    # Subset integrated AnnData per sample id
    ann_adata = adata[adata.obs['sample_id'] == sample_id]
    num_total_cells = ann_adata.shape[0]
    
    # Get condition
    condition = ann_adata.obs['condition'].tolist()[0]
    
    # Open RNA counts and filter by overlaped cells
    input_path = '../data/{0}/outs/filtered_feature_bc_matrix/'.format(sample_id)
    raw_adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', cache=True)
    cell_ids = ['-'.join(c_id.split('-')[:-1]) for c_id in ann_adata.obs.index]
    raw_adata = raw_adata[cell_ids, ann_adata.raw.var.index]
    
    # Get cell type names
    cell_types = np.unique(ann_adata.obs['cell_type'])
    for cell_type in cell_types:
        # Subset annotated data by cell type
        cell_adata = ann_adata[ann_adata.obs['cell_type'] == cell_type]
        # Subset also raw AnnData
        cell_ids = ['-'.join(c_id.split('-')[:-1]) for c_id in cell_adata.obs.index]
        craw_adata = raw_adata[cell_ids]
        
        # Create pseudo bulk profile by summing all counts
        psdo_counts = np.sum(craw_adata.X.toarray(), axis=0)
        
        # Store results
        X.append(psdo_counts)
        col_sample_id.append(sample_id)
        col_condition.append(condition)
        col_cell_type.append(cell_type)
        num_cells = cell_adata.shape[0]
        prop_cells = num_cells / num_total_cells
        col_cell_num.append(num_cells)
        col_cell_prop.append(prop_cells)
        
# Crete AnnData object
pb_adata = AnnData(np.array(X))
pb_adata.obs.index = [sample_id+'_'+cell_type for sample_id,cell_type in zip(col_sample_id, col_cell_type)]
pb_adata.obs['sample_id'] = col_sample_id
pb_adata.obs['condition'] = col_condition
pb_adata.obs['cell_type'] = col_cell_type
pb_adata.obs['cell_num'] = col_cell_num
pb_adata.obs['cell_prop'] = col_cell_prop
pb_adata.var.index = ann_adata.raw.var.index

# Set low counts to 0
pb_adata.X[pb_adata.X <= 2] = 0

def vsn_normalize(arr):
    import logging
    import rpy2.rinterface_lib.callbacks
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import rpy2.robjects as robjects
    
    robjects.globalenv['arr'] = arr
    arr = robjects.r('''
            library(vsn)
            arr <- t(arr)
            arr[is.nan(arr)] = NA
            fit <- vsnMatrix(arr)
            arr <- t(vsn::predict(fit,arr))
            arr[is.na(arr)] <- 0
            arr
            ''')
    return arr

# VSN normalize
for ctype in np.unique(pb_adata.obs['cell_type']):
    # Get expression for ctype
    X = pb_adata.X[pb_adata.obs['cell_type'] == ctype].copy()
    Xdim = X.shape
    # Set zeros to nan
    X[X == 0] = np.nan
    # Filter genes not expressed in almost all samples
    msk = np.sum(np.isnan(X), axis=0) < 2
    X = X[:, msk]
    # VSN normalize
    X_norm = vsn_normalize(X)
    # Set filtered genes expr to 0
    X = np.zeros(Xdim)
    X[:,msk] = X_norm
    pb_adata.X[pb_adata.obs['cell_type'] == ctype] = X
    
# Write to file
pb_adata.write('../qc_data/pseudobulk.h5ad')