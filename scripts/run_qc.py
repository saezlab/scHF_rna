import scanpy as sc
import scanpy.external as sce

import numpy as np
import pandas as pd

import pickle
import os

"""
Script to run basic QC filtering. Stores cell meta data for future plotting and writes new AnnData objects.
"""

meta = {
    'healthy' : ["CK114","CK115","CK139","CK140"],
    'acidosis' : ["CK128"],
    'hf' : ["CK127","CK129","CK135","CK137","CK141"],
    'hf_ckd' : ["CK116","CK126","CK136","CK138"]
}

for condition, samples in meta.items():
    for sample in samples:
        # Read raw data
        print(sample)
        input_path = '../data/{0}/outs/filtered_feature_bc_matrix/'.format(sample)
        adata = sc.read_10x_mtx(input_path,var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        
        # Basic filtering
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # Compute QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        # Compute doublets score
        sce.pp.scrublet(adata, expected_doublet_rate=0.06, verbose=False)
        
        # Set filter values (can be changed)
        mt_thr = 0.5
        gene_qnt = 0.9975
        doublet_thr = adata.uns['scrublet']['threshold']
        
        # Save cell meta data
        df = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score']]
        plot_data = {'mt_thr' : mt_thr,
                     'gene_qnt' : gene_qnt,
                     'doublet_thr' : doublet_thr,
                     'df' : df
                    }
        
        # Filter
        gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)
        msk = (adata.obs.n_genes_by_counts < gene_thr) & \
              (adata.obs.pct_counts_mt < mt_thr) &  \
              (adata.obs.doublet_score < doublet_thr)
        adata = adata[msk, :]
        
        # Save results
        os.makedirs('../plot_data/qc/', exist_ok=True)
        pickle.dump(plot_data, open('../plot_data/qc/{0}.pkl'.format(sample), "wb"))
        
        os.makedirs('../qc_data/', exist_ok=True)
        adata.write('../qc_data/{0}.h5ad'.format(sample))
