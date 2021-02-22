import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd


'''
Open all samples QC processed files, concatenate them and run integration using 
BKNN
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
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, 
                            min_disp=0.5, batch_key='batch')

# Save raw gene expression
adata.raw = adata
adata.uns['hvg']['ngene'] = adata.shape[1]

# Filter HVG. 
# Select top 3000 HVG in as much batches as possible
# Genes need to be HV in at least 2 batches
num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata = adata[:, hvg]
print("Shape:", adata.shape)

# Update QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Compute PCA
sc.tl.pca(adata, svd_solver='arpack')

#sce.pp.harmony_integrate(adata, 'batch', 
#                         adjusted_basis='X_pca', 
#                         max_iter_harmony=30)
# Compute NN
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Integrate using BKNN
sce.pp.bbknn(adata, batch_key='batch')

# Run UMAP
sc.tl.umap(adata)

# Write to file
adata.write('../qc_data/integrated.h5ad')
