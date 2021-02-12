import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd

'''
Open all samples QC processed files, concatenate them and run integration using 
harmony
'''

meta = {
    'healthy' : ["CK114","CK115","CK139","CK140"],
    'acidosis' : ["CK128"],
    'hf' : ["CK127","CK129","CK135","CK137","CK141"],
    'hf_ckd' : ["CK116","CK126","CK136","CK138"]
}

# Open and concatenate all samples
for i, items in enumerate(meta.items()):
    condition, samples = items
    for sample in samples:
        print(sample)
        input_path = '../qc_data/{0}.h5ad'.format(sample)
        # If first sample initalize AnnData object
        if i == 0:
            adata = sc.read_h5ad(input_path)
            adata.obs['sample'] = [sample]*adata.n_obs
            adata.obs['condition'] = [condition]*adata.n_obs
            
        # Else concatenate
        else:
            adata_tmp = sc.read_h5ad(input_path)
            adata_tmp.obs['sample'] = [sample]*adata_tmp.n_obs
            adata_tmp.obs['condition'] = [condition]*adata_tmp.n_obs
            adata = adata.concatenate(adata_tmp)
            
# Identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save raw gene expression
adata.raw = adata

# Filter HVG
adata = adata[:, adata.var.highly_variable]

# Compute PCA
sc.tl.pca(adata, svd_solver='arpack')

# Integrate using harmony
sce.pp.harmony_integrate(adata, 'sample', adjusted_basis='X_pca')

# Compute NN
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Run UMAP
sc.tl.umap(adata)

# Write to file
adata.write('../qc_data/integrated.h5ad')
