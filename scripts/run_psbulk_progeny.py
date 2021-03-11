import scanpy as sc
import numpy as np
import pandas as pd

import progeny
import os

# Read AnnData object
input_path = '../qc_data/pseudobulk.h5ad'
adata = sc.read_h5ad(input_path)

cell_types = np.unique(adata.obs['cell_type'])
conditions = np.unique(adata.obs['condition'])

def lm(data, formula='y ~ x'):
    import statsmodels.api as sm
    
    model = sm.OLS.from_formula(formula, data).fit()
    coeff = model.params[1]
    pval = model.pvalues[1]
    return coeff, pval

def rank_func_feature(adata, name, cond_a, cond_b='healthy', func='progeny'):
    # Filter by conditions to test
    msk = (adata.obs['condition'] == cond_a) | (adata.obs['condition'] == cond_b)
    conds = adata.obs['condition'][msk]
    
    # Reorder conditions
    conds = conds.astype(str).astype('category')
    conds = conds.cat.reorder_categories([cond_a, cond_b])

    # Create input OLS: predict func value (y) by condition (x)
    data = pd.DataFrame()
    data['x'] = conds
    data['y'] = adata.obsm[func][name]

    # Build lm model
    coeff, pval = lm(data)
    
    return coeff, pval


func = 'progeny'
results = []
for cell_type in cell_types:
    subadata = adata[adata.obs['cell_type'] == cell_type]
    progeny.run(subadata, scale=True, organism='Human', top=100)
    names = subadata.obsm[func].columns
    for i,cond_a in enumerate(conditions):
        for j,cond_b in enumerate(conditions):
            if cond_a != cond_b:
                for name in names:
                    coeff, pval = rank_func_feature(subadata, name, cond_a, cond_b, func=func)
                    result = [cell_type, cond_a, cond_b, name, coeff, pval]
                    results.append(result)

# Transform to df
results = pd.DataFrame(results, columns=['cell_type', 'ref', 'cond', 'name', 'coeff', 'pvals']).sort_values('pvals')

# Save
os.makedirs('../plot_data/func/', exist_ok=True)
results.to_csv('../plot_data/func/progeny.csv', index=False)