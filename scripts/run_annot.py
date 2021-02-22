import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

'''Annotate the obtianed clusters by a specified resolution'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Select desired resolution
res = 'leiden_res_1.0'

# Write annotated clusters
annot_clus = {
    0 : 'cardiomyocyte',
    1 : 'endothelial',
    2 : 'fibroblast',
    3 : 'fibroblast',
    4 : 'monocytes',
    5 : 'pericyte',
    6 : 'fibroblast',
    7 : 'T-cells',
    8 : 'monocytes',
    9 : 'endothelial',
    10 : 'cardiomyocyte',
    11 : 'vSMCs',
    12 : 'unknown',
    13 : 'mast_cells',
    14 : 'adipocytes',
    15 : 'neuronal',
    16 : 'lymphatic_endo',
    17 : 'unknown',
    18 : 'unknown',
    19 : 'endothelial'
}

# Add metada-data
adata.obs['cell_type'] = [annot_clus[int(clus)] for clus in adata.obs[res]]

# Write to file
adata.write('../qc_data/integrated.h5ad')