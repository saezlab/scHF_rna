import scanpy as sc
import numpy as np
import pandas as pd

from utils import cos_dis, get_group_sum_dis
from scipy import stats
from statsmodels.stats.multitest import multipletests

def get_cosine_distances(input_path, sign_path, n_features=500):
    # Read inputs
    adata = sc.read_h5ad(input_path)
    sign = pd.read_csv(sign_path)

    # Get unique conditions, cell_types, contrasts
    conditions = np.unique(adata.obs['condition'])
    cell_types = np.unique(adata.obs['cell_type'])
    contrasts = np.unique(sign['contrast'])

    # Compute angular distances for each condition
    df = []
    for contrast in contrasts:
        cond_a, cond_b = contrast.split('vs')
        col_cell_types = []
        col_dists = []
        for cell_type in cell_types:

            # Get a subset per contrast and cell type
            adata_a = adata[(adata.obs['condition']==cond_a)&(adata.obs['cell_type']==cell_type)]
            adata_b = adata[(adata.obs['condition']==cond_b)&(adata.obs['cell_type']==cell_type)]

            # Get union of top n significant features per contrasts for a given cell type
            ftrs = sign[sign['cell_type']==cell_type]
            ftrs_lst = ftrs[ftrs['contrast']==contrast].sort_values('pvals').head(n_features)['names'].values

            ftrs_lst = np.unique(ftrs_lst)
            print(cell_type, len(ftrs_lst))
            if len(ftrs_lst) != 0:
                adata_a = adata_a[:,ftrs_lst]
                adata_b = adata_b[:,ftrs_lst]

            #  Compute mean within groups as norm factor
            sum_a, ndists_a = get_group_sum_dis(adata_a)
            sum_b, ndists_b = get_group_sum_dis(adata_b)
            norm = (sum_a + sum_b) / (ndists_a + ndists_b)
            if norm == 0.0:
                norm = 1

            # Compute all possible distances
            for row_a in adata_a.X.toarray():
                    for row_b in adata_b.X.toarray():
                        # Compute angular distance
                        if len(ftrs_lst) != 0:
                            dis = cos_dis(row_a,row_b)/norm
                        else:
                            dis = np.nan
                        row = [contrast, cell_type, dis]
                        df.append(row)
    df = pd.DataFrame(df, columns=['contrast', 'cell_type', 'dist'])
    return df

# RNA distances
input_path = '../qc_data/pseudobulk.h5ad'
sign_path = '../plot_data/deg/psbulk_deg.csv'
print('RNA')
rna_df = get_cosine_distances(input_path, sign_path, n_features=500)
rna_df['type'] = 'rna'

# ATAC distances
input_path = '../qc_data/atac_peaks.h5ad'
sign_path = '../plot_data/deg/psbulk_dap.csv'
print('ATAC')
atac_df = get_cosine_distances(input_path, sign_path, n_features=500)
atac_df['type'] = 'atac'

# Merge
df = pd.concat([rna_df, atac_df])

# Compute significance contrasts
table = []
for contrast in np.unique(df['contrast']):
    dist = df[df['contrast'] == contrast].dist.values
    dist = dist[~np.isnan(dist)]
    _, pval = stats.ttest_1samp(dist, 1, alternative='greater')
    table.append([contrast, pval])
table = pd.DataFrame(table, columns=['contrast', 'pval'])
    
# Correct by FDR
_, pvals_adj, _, _ = multipletests(table['pval'].values, alpha=0.05, method='fdr_bh')
table['adj_pval'] = pvals_adj
table = table.sort_values('pval')

# Save
table.to_csv('../tables/shifts_contrasts.csv', index=False)

# Compute significance cell types
table = []
for cell_type in np.unique(df['cell_type']):
    for contrast in np.unique(df['contrast']):
        dist = df[(df['cell_type'] == cell_type)&(df['contrast'] == contrast)].dist.values
        dist = dist[~np.isnan(dist)]
        _, pval = stats.ttest_1samp(dist, 1, alternative='greater')
        table.append([cell_type, contrast, pval])
table = pd.DataFrame(table, columns=['cell_type', 'contrast', 'pval'])

# Correct by FDR
_, pvals_adj, _, _ = multipletests(table['pval'].values, alpha=0.05, method='fdr_bh')
table['adj_pval'] = pvals_adj
table = table.sort_values('pval')

# Save
table.to_csv('../tables/shifts_celltypes.csv', index=False)

# Save
df.to_csv('../plot_data/var_shifts/shifts.csv', index=False)
