import scanpy as sc
import numpy as np
import pandas as pd


healthy = ["CK114","CK115","CK139","CK140"]
acidosis = ["CK128"]
hf = ["CK127","CK129","CK135","CK137","CK141"]
hf_ckd = ["CK116","CK126","CK136","CK138"]

samples = healthy + acidosis + hf + hf_ckd
qc_stats = []

for sample in samples[:2]:
    print(sample)

    # Read sample
    input_path = '../../../bq_shared/scellHF/scheart_data_red/{0}/outs/filtered_feature_bc_matrix/'.format(sample)
    adata = sc.read_10x_mtx(
        input_path,  
        var_names='gene_symbols',
        cache=False)
    adata.var_names_make_unique()
    
    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Store metadata before filtering
    before_obs = adata.obs.copy()

    # Filter
    n_genes_by_counts_filter = np.quantile(adata.obs.n_genes_by_counts, 1-0.005)
    mt_per_filter = 5
    adata = adata[adata.obs.n_genes_by_counts < n_genes_by_counts_filter, :]
    adata = adata[adata.obs.pct_counts_mt < mt_per_filter, :]
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    after = adata.obs.copy()
