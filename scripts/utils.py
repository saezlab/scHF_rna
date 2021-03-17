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