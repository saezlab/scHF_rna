import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

'''Annotate the obtianed clusters by a specified resolution'''

# Read AnnData object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Select desired resolution
res = 'leiden_res_1.25'

# Write annotated clusters
annot_clus = {
    0 : 'endothelial',
    1 : 'cardiomyocyte',
    2 : 'macrophages',
    3 : 'pericyte',
    4 : 'fibroblast',
    5 : 'fibroblast',
    6 : 'fibroblast',
    7 : 'fibroblast',
    8 : 'cardiomyocyte',
    9 : 'macrophages',
    10 : 'endothelial',
    11 : 'cardiomyocyte',
    12 : 'T-cells',
    13 : 'cardiomyocyte',
    14 : 'cardiomyocyte',
    15 : 'vSMCs',
    16 : 'endothelial',
    17 : 'mast_cells',
    18 : 'T-cells',
    19 : 'adipocytes',
    20 : 'macrophages',
    21 : 'neuronal',
    22 : 'endothelial',
    23 : 'lymphatic_endo',
    24 : 'fibroblast',
    25 : 'unknown',
    26 : 'unknown',
    27 : 'fibroblast'
}

# Add metada-data
adata.obs['cell_type'] = [annot_clus[int(clus)] for clus in adata.obs[res]]

# Write to file
adata.write('../qc_data/integrated.h5ad')