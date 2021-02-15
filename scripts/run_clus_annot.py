import scanpy as sc
import numpy as np
import pandas as pd

from SCCAF import SCCAF_assessment, plot_roc, SCCAF_optimize_all

''' '''


# Read integrated object
input_path = '../qc_data/integrated.h5ad'
adata = sc.read_h5ad(input_path)

# Overcluster
sc.tl.leiden(adata, resolution=1.5, key_added='L2_Round0')
SCCAF_optimize_all(ad=adata, plot=False, min_acc=0.9, prefix = 'L2', use='umap')