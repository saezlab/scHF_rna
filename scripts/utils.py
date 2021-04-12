import numpy as np
import pandas as pd

def vsn_normalize(arr):
    '''
    Normalizes expression array (samples x genes) by variance stabilization (vsn).
    '''
    import logging
    import rpy2.rinterface_lib.callbacks
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import rpy2.robjects as robjects
    
    robjects.globalenv['arr'] = arr
    arr = robjects.r('''
            library(vsn)
            arr <- t(arr)
            arr[is.nan(arr)] = NA
            fit <- vsnMatrix(arr)
            arr <- t(vsn::predict(fit,arr))
            arr[is.na(arr)] <- 0
            arr
            ''')
    return arr

def cos_sim(a,b):
    return np.min([np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)),1])

def ang_dis(a,b):
    return np.arccos(cos_sim(a,b)) / np.pi

def get_group_sum_dis(sub_adata):
    # Compute all cumulative distances in an AnnData object
    cum_dis = 0
    for i,row_a in enumerate(sub_adata.X.toarray()):
        for j, row_b in enumerate(sub_adata.X.toarray()):
            if j < i:
                # Compute angular distance
                dis = ang_dis(row_a,row_b)
                cum_dis += dis
    return cum_dis

def lm(data, formula='y ~ x'):
    '''
    Builds a one covariable linear model and returns its coeff and pval
    '''
    import statsmodels.api as sm
    model = sm.OLS.from_formula(formula, data).fit()
    coeff = model.params[1]
    pval = model.pvalues[1]
    return coeff, pval

def rank_func_feature(adata, name, cond_a, cond_b, func):
    '''
    Computes differences of functional features across two conditions
    '''
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

def get_db(path='../data/canonical_pathways.gmt'):
    # Get descriptin for each term id
    descr_df = pd.read_csv('../data/msigdbr_id_to_descr.csv')
    descr_dict = {k:v for k,v in zip(descr_df['gs_name'], descr_df['gs_description'])}
    db = []
    # Read GMT file
    with open(path, "r") as f:
        for line in f.readlines():
            line = line.rstrip().split('\t')
            g_set_id = line[0]
            descr = descr_dict[g_set_id]
            g_set = line[2:]
            db.append([g_set_id, descr, g_set])
    return db

def ORA(g_lst, db):
    from statsmodels.stats.multitest import multipletests
    import scipy.stats as stats
    
    # Select genes that are in the given db
    g_lst = set(g_lst)
    g_db = set()
    for row in db:
        for gname in row[2]:
            g_db.add(gname)
    g_lst = g_lst & g_db
    
    num_g_lst = len(g_lst)
    num_g_db = len(g_db)
    
    col_g_set_id = []
    col_descr = []
    col_num_set = []
    col_num_set_lst = []
    col_gnames = []
    col_pvalues = []
    
    for row in db:
        g_set_id, descr, g_set = row
        g_set = set(g_set)
        num_set = len(g_set)
        num_set_lst = len(g_set & g_lst)
        
        col_g_set_id.append(g_set_id)
        col_descr.append(descr)
        col_num_set.append(num_set)
        col_num_set_lst.append(num_set_lst)
        col_gnames.append(list(g_set))
        
        pval = stats.hypergeom.sf(k=num_set_lst-1, M=num_g_db, n=num_set, N=num_g_lst, loc=0)
        col_pvalues.append(pval)
    
    df = pd.DataFrame()
    df['g_set_id'] = col_g_set_id
    df['descr'] = col_descr
    df['num_set'] = col_num_set
    df['num_set_lst'] = col_num_set_lst
    df['gnames'] = col_gnames
    df['pvalues'] = col_pvalues
    df['adj_pvalues'] = multipletests(col_pvalues, method='fdr_bh')[1]
        
    return df

def read_SummarizedExperiment(file_path, atlas_name, ann_cells):
    '''
    Gets expression df from SUmmarizedExperiment list
    '''
    import logging
    import rpy2.rinterface_lib.callbacks
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import rpy2.robjects as robjects
    
    import scanpy as sc
    from anndata import AnnData
    
    robjects.globalenv['atlas_name'] = atlas_name
    robjects.globalenv['ann_cells'] = ann_cells
    robjects.globalenv['file_path'] = file_path
    expr_df = robjects.r('''
            data <- readRDS('../qc_data/hca_pseudobulk.rds')
            expr_df <- data[[atlas_name]][[ann_cells]]@assays@data@listData[["sum"]]
            as.data.frame(expr_df)
            ''')
    
    # Transpose and as int
    expr_df = expr_df.astype(np.int).T
    
    # Generate AnnData
    atlas = AnnData(np.array(expr_df), obs=pd.DataFrame(index=expr_df.index), 
                    var=pd.DataFrame(index=expr_df.columns))
    
    # Filter genes
    sc.pp.filter_genes(atlas, min_cells=1)
    
    return atlas

def get_design(df, col):
    x = df.pivot(columns='condition', values='condition').fillna(0)
    x[x != 0] = 1
    return x

def get_contrast(design, contr_dict):
    df = pd.DataFrame(columns=contr_dict.keys(), index=design.columns)
    for contrast in contr_dict:
        cont, ref = contr_dict[contrast]
        df.loc[cont,contrast] = 1
        df.loc[ref,contrast] = -1
    return df.fillna(0)

def limma_fit(X, design, contr_matrix):
    '''
    Fits limma model
    '''
    import logging
    import rpy2.rinterface_lib.callbacks
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import rpy2.robjects as robjects
    
    robjects.globalenv['X'] = X
    robjects.globalenv['design'] = design
    robjects.globalenv['contr_matrix'] = contr_matrix
    x = robjects.r('''
            library(limma)
            fit <- lmFit(X, design)
            fit2 <- contrasts.fit(fit, contr_matrix)
            fit2 <- eBayes(fit2)
            coefs <- as.data.frame(fit2$coefficients)
            pvals <- as.data.frame(fit2$p.value)
            list(coefs, pvals)
            ''')
    coefs, pvals = x[0], x[1]
    df = coefs.melt(ignore_index=False, value_name='logfoldchanges', var_name='contrast')
    df['pvals'] = pvals.melt(ignore_index=False, value_name='pvals', var_name='contrast')['pvals']
    df = df.reset_index()
    df.columns = ['names', 'contrast', 'logfoldchanges', 'pvals']
    
    return df