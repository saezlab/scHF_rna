import scanpy as sc
import numpy as np
import pandas as pd

import os
import pickle

from utils import cos_sim, ang_dis, get_group_sum_dis

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

# Load DEG
df_deg = pd.read_csv('../plot_data/deg/psbulk_deg.csv')

# Get unique conditions and cell_types
conditions = np.unique(adata.obs['condition'])
cell_types = np.unique(adata.obs['cell_type'])

plot_data = dict()

for i,cond_a in enumerate(conditions):
    for j,cond_b in enumerate(conditions):
        if i < j:
            col_cell_types = []
            col_dists = []
            for cell_type in cell_types:
                # Get a subset per condition and cell type
                adata_a = adata[(adata.obs['condition']==cond_a)&(adata.obs['cell_type']==cell_type)]
                adata_b = adata[(adata.obs['condition']==cond_b)&(adata.obs['cell_type']==cell_type)]
                
                # Get union of top 100 genes that are sign DEG per contrasts for a given cell type
                deg = df_deg[df_deg['cell_type']==cell_type]
                deg = np.unique([deg[deg['contrast']==contrast].sort_values('pvals').head(100)['names'].tolist() \
                                 for contrast in np.unique(deg['contrast'])])
                print(cell_type, len(deg))
                adata_a = adata_a[:,deg]
                adata_b = adata_b[:,deg]
                
                #  Compute mean within groups as norm factor
                norm = (get_group_sum_dis(adata_a) + get_group_sum_dis(adata_b)) / \
                (adata_a.shape[0] + adata_b.shape[0])
                
                if norm == 0.0:
                    norm = 1
                # Compute all possible interactions
                for row_a in adata_a.X.toarray():
                        for row_b in adata_b.X.toarray():
                            # Compute angular distance
                            dis = ang_dis(row_a,row_b)/norm
                            col_dists.append(dis)
                            col_cell_types.append(cell_type)
            # Save to dict
            df = pd.DataFrame()
            df['cell_type'] = col_cell_types
            df['dists'] = col_dists
            plot_data['{0}.vs.{1}'.format(cond_a,cond_b)] = df

# Save to file
os.makedirs('../plot_data/angdist/', exist_ok=True)
pickle.dump(plot_data, open('../plot_data/angdist/plot_data.pkl', "wb"))
