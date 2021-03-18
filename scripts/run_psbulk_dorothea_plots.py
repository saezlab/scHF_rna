import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from plotting import dotplot


# Read tested TFs
results = pd.read_csv('../plot_data/func/dorothea.csv')
results['ctype.cond'] = results['cell_type'] + ' ' + results['cond']
refs = np.unique(results['ref'])

# Filter by cell type
cell_type = 'fibroblast'

for ref in refs:
    # Filter by a reference condition
    df = results[(results['ref'] == ref)]
    
    # Get significant elements
    sign = df[df['pvals'] < 0.05]
    sign = sign[sign['cell_type'] == cell_type].head(50)
    names = np.unique(sign['name'].tolist())
    ctconds = np.unique(sign['ctype.cond'].tolist())
    
    # Filter by combination of significant cell type and condition
    msk_ctconds = np.array([ctcond in ctconds for ctcond in df['ctype.cond'].tolist()])
    msk_names = np.array([name in names for name in df['name'].tolist()])
    df = df[msk_ctconds * msk_names]
    msk = df['pvals'] > 0.05
    
    # Hide non significant
    df['pvals'].loc[msk] = 1
    df['coeff'].loc[msk] = 0
    
    # Plot dot plot
    fig = dotplot(ref, df)
    fig.savefig('../plots/psbulk_dorothea_{0}.png'.format(ref), bbox_inches='tight')