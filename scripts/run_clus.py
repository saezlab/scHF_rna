import scanpy as sc
import numpy as np
import pandas as pd

'''Clusters cells in different leiden resolutions'''

# Read integrated object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

scores = []
# Sequence of clustering resolutions
for res in np.arange(0.5, 1.05, 0.1):
    res = np.around(res, 2)
    res_name = 'leiden_res_{0}'.format(res)
    print(res_name)
    # Cluster based on a specific resolution
    sc.tl.leiden(adata, resolution=res, key_added=res_name)


# Write to file
adata.write('../qc_data/integrated.h5ad')