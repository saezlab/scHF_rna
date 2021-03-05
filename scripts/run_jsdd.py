import scanpy as sc
import numpy as np
import pandas as pd

import pickle
import os

from sklearn.manifold import MDS

from scipy.spatial.distance import jensenshannon



# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

# Load DEG
df_deg = pd.read_csv('../plot_data/deg/deg.csv')

# Get list samples and cell types
samples_ids = np.unique(adata.obs['sample_id'])
conditions = [adata.obs.condition[adata.obs['sample_id'] == sample_id].tolist()[0] for sample_id in samples_ids]
cell_types = np.unique(adata.obs['cell_type'])


t_jsd_ds = []
t_props = []
t_coords = []

for cell_type in cell_types:
    # Subset by cell type
    sub_adata = adata[adata.obs['cell_type'] == cell_type]
    
    # Init JSDD matrix and props
    jsd_ds = []
    props = []
    
    # Get genes that are sign DEG in all conditions for this cell type
    deg = df_deg[df_deg['cell_type']==cell_type]
    deg = np.unique([deg[deg['condition']==cond].sort_values('pvals').head(100)['names'].tolist() \
                     for cond in np.unique(deg['condition'])])
    print(cell_type, len(deg))
    
    # Compute JDDS for all combinations
    for sample_a in samples_ids:
        for sample_b in samples_ids:
            msk_a = np.array(sub_adata.obs['sample_id'] == sample_a)
            msk_b = np.array(sub_adata.obs['sample_id'] == sample_b)
            
            # If one of the samples has no cells skip
            if np.sum(msk_a) == 0 or np.sum(msk_b) == 0:
                jsd_d = np.nan
                w = 0
            else:
                # Get sample objects
                x_a = sub_adata[msk_a]
                x_b = sub_adata[msk_b]
                
               # Filter by DEG
                x_a = x_a[:,deg].X.toarray()[0]
                x_b = x_b[:,deg].X.toarray()[0]
                
                # Get min proportion
                w_a = sub_adata.obs['cell_prop'][msk_a]
                w_b = sub_adata.obs['cell_prop'][msk_b]
                w = np.min([w_a, w_b])
                
                # Compute JSDD and normalize by min prop
                jsd_d = jensenshannon(x_a, x_b)
                
            # Store result
            jsd_ds.append(jsd_d)
            props.append(w)
    
    # Reshape to matrix and store
    jsd_ds = np.array(jsd_ds).reshape(len(samples_ids), len(samples_ids))
    props = np.array(props).reshape(len(samples_ids), len(samples_ids))
    t_jsd_ds.append(jsd_ds)
    t_props.append(props)
    
    # Compute MDS coordinates of non nans and store
    msk = ~np.isnan(np.diag(jsd_ds))
    coords = np.full((len(samples_ids),2), np.nan)
    coords[msk] = MDS(dissimilarity='precomputed', random_state=0).fit_transform(jsd_ds[msk,][:,msk])    
    t_coords.append(coords)

# Compute final JSDD for all cell types
all_jsd_ds = np.sum(np.nan_to_num(t_jsd_ds)*t_props, axis=0)
all_coords = MDS(dissimilarity='precomputed', random_state=0).fit_transform(all_jsd_ds)

# Store results in a dict
plot_data = {
    'samples_ids' : np.array(samples_ids),
    'conditions' : np.array(conditions),
    'cell_types' : np.hstack(['All', cell_types]),
    'jsds' : [all_jsd_ds] + t_jsd_ds,
    'coords' : [all_coords] + t_coords
}

# Save results
os.makedirs('../plot_data/jsd/', exist_ok=True)
pickle.dump(plot_data, open('../plot_data/jsd/plot_data.pkl', "wb"))