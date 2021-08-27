import scanpy as sc
import numpy as np
import pandas as pd

import os
from anndata import AnnData

path = '../data/ATAC/DiffPeaksBulkMatrix'
files = os.listdir(path)

dap = []
for file in files:
    # Get cell type and contrast from file name
    contrast = file.split('_')
    if contrast[-1] == 'a.csv':
        cond2 = 'HF-A'
    elif contrast[-1] == 'ckd.csv':
        cond2 = 'HF-CKD'
    if contrast[0] == 'mast':
        cell_type = 'mast_cells'
        cond1 = contrast[2]
        if cond1 == 'control':
            cond1 = 'Healthy'
    else:
        cell_type = contrast[0]
        cond1 = contrast[1]
        if cond1 == 'control':
            cond1 = 'Healthy'
        else:
            if contrast[2] == 'a':
                cond1 = 'HF-A'
            elif contrast[2] == 'ckd':
                cond1 = 'HF-CKD'

    if cond1 == 'hf':
        cond1 = 'HF-A'
    contrast = cond2 + '-' + cond1
    names = pd.read_csv(os.path.join(path, file))
    #names = np.unique(names.peak_pos)
    df = pd.DataFrame(columns=['names', 'contrast', 'cell_type'])
    df['names'] = names.peak_pos
    df['contrast'] = contrast
    df['cell_type'] = cell_type
    df['pvals'] = names.PValue
    dap.append(df)
df_dap = pd.concat(dap)

df_dap.to_csv('../plot_data/deg/psbulk_dap.csv', index=False)