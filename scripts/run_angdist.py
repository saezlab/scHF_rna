import scanpy as sc
import numpy as np
import pandas as pd

import os
import pickle


def cos_sim(a,b):
    return np.min([np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)),1])

def ang_dis(a,b):
    return np.arccos(cos_sim(a,b)) / np.pi

def get_group_sum_dis(sub_adata):
    # Compute all cumulative distances in an AnnData object
    cum_dis = 0
    for i,row_a in enumerate(sub_adata.X.toarray()):
        for j, row_b in enumerate(sub_adata.X.toarray()):
            if j < i:
                # Compute angular distance
                dis = ang_dis(row_a,row_b)
                cum_dis += dis
    return cum_dis

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

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