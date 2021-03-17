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
    0 : 'macrophages',
    1 : 'cardiomyocyte',
    2 : 'fibroblast',
    3 : 'endothelial',
    4 : 'pericyte',
    5 : 'fibroblast',
    6 : 'fibroblast',
    7 : 'T-cells',
    8 : 'fibroblast',
    9 : 'cardiomyocyte',
    10 : 'cardiomyocyte',
    11 : 'endothelial',
    12 : 'vSMCs',
    13 : 'endothelial',
    14 : 'mast_cells',
    15 : 'adipocytes',
    16 : 'neuronal',
    17 : 'endothelial',
    18 : 'macrophages',
    19 : 'lymphatic_endo',
    20 : 'unknown',
    21 : 'unknown',
    22 : 'unknown'
}

# Add metada-data
adata.obs['cell_type'] = [annot_clus[int(clus)] for clus in adata.obs[res]]

# Write to file
adata.write('../qc_data/integrated.h5ad')