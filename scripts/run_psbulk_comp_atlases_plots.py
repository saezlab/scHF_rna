# Compare cell types agains other atlases

# Must first generate a unique profile per cell type

import scanpy as sc
import numpy as np
import pandas as pd

from anndata import AnnData
from scipy.spatial.distance import jensenshannon
from utils import read_SummarizedExperiment

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

# Get unique cell types
cell_types = np.unique(adata.obs['cell_type'])

X = []
col_cell_type = []

for cell_type in cell_types:
    # Subset by cell type
    subadata = adata[adata.obs['cell_type'] == cell_type]
    
    # Get sum of all genes
    psbulk_counts = np.sum(subadata.raw.X, axis=0)
    
    # Append results
    X.append(psbulk_counts)
    col_cell_type.append(cell_type)

X = np.array(X)
adata = AnnData(X, var=pd.DataFrame(index=adata.var.index))
adata.obs.index = cell_types
adata.obs['cell_type'] = col_cell_type

# Store raw counts
adata.raw = adata


# Filter genes
sc.pp.filter_genes(adata, min_cells=1)

# Read atlas
atlas_name = 'cells'
ann_cells = 'cell_states'

for atlas_name in ['all', 'cells', 'nuclei']:
    for ann_cells in ['cell_type', 'cell_states']:

        atlas = read_SummarizedExperiment('../qc_data/hca_pseudobulk.rds', atlas_name=atlas_name, ann_cells=ann_cells)

        # Subset genes
        msk = atlas.var.index.intersection(adata.var.index)
        fadata = adata[:,msk]
        atlas = atlas[:,msk]

        # Normalize and trasnform per cell type
        sc.pp.normalize_total(fadata, target_sum=1e4)
        sc.pp.log1p(fadata)

        sc.pp.normalize_total(atlas, target_sum=1e4)
        sc.pp.log1p(atlas)

        # JSD distances
        X = []
        for ctype_a in fadata.X:
            row = []
            for ctype_b in atlas.X:
                row.append(jensenshannon(ctype_a, ctype_b))
            X.append(row)
        X = np.array(X)
        
        ######
        min_cols_row_pos = np.argmin(X, axis=0)
        min_cols_order = np.argsort(min_cols_row_pos)
        X = X[:,min_cols_order]
        xlabels = atlas.obs.index[min_cols_order]
        ######
        
        import seaborn as sns
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1,1, figsize=(16, 6), dpi=150)
        sns.heatmap(X, 
                    cmap='Blues_r', 
                    xticklabels=xlabels, 
                    yticklabels=adata.obs.index, 
                    cbar_kws={"shrink": .75, 'label': 'JSD distances'},
                    robust=True,
                    square=True,
                    ax=ax
                    )

        ax.set_title('Human Heart Atlas ({0} | {1}) vs Heart Failure Atlas'.format(atlas_name, ann_cells), fontsize=16)
        fig.set_facecolor('white')
        fig.savefig('../plots/atlas_{0}_{1}'.format(atlas_name, ann_cells), bbox_inches='tight')
