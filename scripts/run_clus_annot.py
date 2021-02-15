import scanpy as sc
import numpy as np
import pandas as pd

from sklearn.metrics import silhouette_score

'''Clusters cells in different leiden resolutions and finds the one that maximises
the average silhouette scores'''

# Read integrated object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

scores = []
# Sequence of clustering resolutions
for res in np.arange(0.1, 1.0, 0.1):
    res = np.around(res, 2)
    res_name = 'leiden_res_{0}'.format(res)
    print(res_name)
    # Cluster based on a specific resolution
    sc.tl.leiden(adata, resolution=res, key_added=res_name)
    # Compute its silhouette score
    score = silhouette_score(adata.obsp['distances'], adata.obs[res_name])
    # Store score
    scores.append([res_name, score])

# Sort by score
scores = np.array(sorted(scores, key=lambda x: x[1], reverse=True))

# Store in AnnData object
adata.uns['silhouette_scores'] = scores
opt_res = scores[-1][0]
adata.obs['opt_clus'] = adata.obs[opt_res]

# Write to file
adata.write('../qc_data/integrated.h5ad')