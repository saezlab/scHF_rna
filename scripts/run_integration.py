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
adatas = []
for condition, samples in meta.items():
    for sample_id in samples:
        print(sample_id)
        input_path = '../qc_data/{0}.h5ad'.format(sample_id)
        adata = sc.read_h5ad(input_path)
        adata.obs['sample_id'] = [sample_id]*adata.n_obs
        adata.obs['condition'] = [condition]*adata.n_obs
        adatas.append(adata)

# Merge objects and delete list
adata = adatas[0].concatenate(adatas[1:])
del adatas

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata, batch_key='batch')

# Save raw gene expression
adata.raw = adata

# Filter HVG. 
# HVG must be in at least the 75% of the batches
num_batches = len(np.unique(adata.obs.batch))
per_batches = 0.75 # Can be changed
min_batch = np.ceil(num_batches * per_batches)
adata.var.highly_variable = adata.var.highly_variable_nbatches > min_batch
adata = adata[:, adata.var.highly_variable]

# Compute PCA
sc.tl.pca(adata, svd_solver='arpack')

# Integrate using harmony
sce.pp.harmony_integrate(adata, 'sample_id', 
                         adjusted_basis='X_pca', 
                         max_iter_harmony=30)

# Compute NN
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Run UMAP
sc.tl.umap(adata)

# Write to file
adata.write('../qc_data/integrated.h5ad')
