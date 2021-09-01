import scanpy as sc
import numpy as np
import pandas as pd

import pickle
import os

from sklearn.manifold import MDS
from scipy.spatial.distance import jensenshannon

from sklearn.metrics import silhouette_score


def run_MDS(adata, df_deg, samples_ids, cell_types, conditions):
    t_jsd_ds = []
    t_props = []
    t_coords = []

    df = []

    for cell_type in cell_types:
        # Subset by cell type
        sub_adata = adata[adata.obs['cell_type'] == cell_type]

        # Init JSDD matrix and props
        jsd_ds = []
        props = []

        # Get genes that are sign DEG in all conditions for this cell type
        deg = df_deg[df_deg['cell_type']==cell_type]
        deg = np.unique(deg.sort_values(['contrast', 'pvals']).groupby('contrast').head(500).names)
        print(cell_type, len(deg))
        sub_adata = sub_adata[:,deg]

        # Compute JDDS for all combinations
        for sample_a in samples_ids:
            for sample_b in samples_ids:
                msk_a = np.array(sub_adata.obs['sample_id'] == sample_a)
                msk_b = np.array(sub_adata.obs['sample_id'] == sample_b)

                # If one of the samples has no cells or genes skip
                if (np.sum(msk_a) == 0 or np.sum(msk_b) == 0) or (sub_adata.shape[1] == 0):
                    jsd_d = np.nan
                    w = 0
                else:
                    # Get sample objects
                    x_a = sub_adata[msk_a]
                    x_b = sub_adata[msk_b]

                    # Filter by DEG
                    x_a = x_a.X.toarray()[0]
                    x_b = x_b.X.toarray()[0]

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
        s_score = np.nan
        if not np.all(np.isnan(jsd_ds)):
            coords[msk] = MDS(dissimilarity='precomputed', random_state=0).fit_transform(jsd_ds[msk,][:,msk])
            s_score = silhouette_score(coords[msk], labels=conditions[msk])
        t_coords.append(coords)

        # Store in df
        for i,sample_id in enumerate(samples_ids):
            cond = conditions[i]
            x, y = coords[i]
            row = [cell_type, cond, sample_id, x, y, s_score]
            df.append(row)

    # Compute final JSDD for all cell types
    all_jsd_ds = np.sum(np.nan_to_num(t_jsd_ds)*t_props, axis=0)
    all_coords = MDS(dissimilarity='precomputed', random_state=0).fit_transform(all_jsd_ds)
    s_score = silhouette_score(all_coords, labels=conditions)
    for i,sample_id in enumerate(samples_ids):
        cond = conditions[i]
        x, y = all_coords[i]
        row = ['All', cond, sample_id, x, y, s_score]
        df.append(row)

    cords = pd.DataFrame(df, columns=['cell_type','condition','sample_id','x','y', 's_score'])
    dists = t_jsd_ds + [all_jsd_ds]
    return cords, np.array(dists)


def center_0_to_1(arr):
    return (arr - np.min(arr)) / (np.max(arr) - np.min(arr))


# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)
df_deg = pd.read_csv('../plot_data/deg/psbulk_deg.csv')

samples_ids = np.unique(adata.obs['sample_id'])
cell_types = np.unique(adata.obs['cell_type'])
conditions = np.array([adata.obs.condition[adata.obs['sample_id'] == sample_id].tolist()[0] for sample_id in samples_ids])


rna_cords, rna_dists = run_MDS(adata, df_deg, samples_ids, cell_types, conditions)
rna_cords['type'] = 'RNA'

# Read AnnData object
input_path = '../qc_data/atac_peaks.h5ad'
adata = sc.read_h5ad(input_path)
df_deg = pd.read_csv('../plot_data/deg/psbulk_dap.csv')

atac_cords, atac_dists = run_MDS(adata, df_deg, samples_ids, cell_types, conditions)
atac_cords['type'] = 'ATAC'

rna_dists[np.isnan(rna_dists)] = 0
atac_dists[np.isnan(atac_dists)] = 0

rna_dists = center_0_to_1(rna_dists)
atac_dists = center_0_to_1(atac_dists)

both_dists = rna_dists + atac_dists

cell_types = np.hstack([cell_types, 'All'])

df = []
for j,arr in enumerate(both_dists):
    cell_type = cell_types[j]
    coords = MDS(dissimilarity='precomputed', random_state=0).fit_transform(arr)
    s_score = silhouette_score(coords, labels=conditions)
    for i,sample_id in enumerate(samples_ids):
            cond = conditions[i]
            x, y = coords[i]
            row = [cell_type, cond, sample_id, x, y, 'BOTH', s_score]
            df.append(row)
            
df = pd.concat([
    rna_cords,
    atac_cords,
    pd.DataFrame(df, columns=['cell_type','condition','sample_id','x','y','type','s_score'])
])

os.makedirs('../plot_data/jsd/', exist_ok=True)
df.to_csv('../plot_data/jsd/jsd.csv', index=False)