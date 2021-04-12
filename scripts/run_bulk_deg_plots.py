import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from plotting import dotplot


# Read DEG df
input_path = '../plot_data/deg/psbulk_deg.csv'
df = pd.read_csv(input_path).set_index('names')
df['logfoldchanges'] = -df['logfoldchanges']
df['contrast'] = ['-'.join(contrast.split('-')[::-1]) for contrast in df['contrast']]

input_path = '../plot_data/deg/reheat_deg.csv'
rh = pd.read_csv(input_path).head(500)
rh.columns = ['names', 'logfoldchanges']
rh = rh.set_index('names')

input_path = '../plot_data/deg/bulk_deg.csv'
bf = pd.read_csv(input_path).groupby('contrast').head(500).set_index('names')

def get_contrast_corr(df_a, df_b, col_name='logfoldchanges'):
    from scipy import stats
    
    contrasts = np.unique(df_a.contrast)
    cell_types = np.unique(df_a.cell_type)
    
    c_label = []
    x_label = []
    y_label = []
    s_label = []
    
    for contrast in contrasts:
        contrast_df = df_a[df_a.contrast == contrast]
        for ctype in cell_types:
            contrast_ctype_df = contrast_df[contrast_df.cell_type == ctype]
            common_genes = contrast_ctype_df.index.intersection(df_b.index)
            t_df_b = df_b.loc[common_genes]
            contrast_ctype_df = contrast_ctype_df.loc[common_genes]
            x, y = contrast_ctype_df[col_name].values, t_df_b[col_name].values
            if len(x) != 0 or len(y) != 0:
                corr, _ = stats.spearmanr(x, y)
            else:
                corr = 0
            x_label.append(ctype)
            y_label.append(contrast)
            s_label.append(len(common_genes))
            c_label.append(corr)
    return np.array(x_label), np.array(y_label), np.array(s_label), np.array(c_label)

# ReHeat
x_label, y_label, s_label, c_label = get_contrast_corr(df, rh)
fig = dotplot(
        title='scHF vs ReHeat', 
        x=x_label, 
        y=y_label, 
        c=c_label, 
        s=s_label, 
        size_title='Overlap genes', 
        color_title='Spearman corr',
        num=0.75,
        figsize=(9,4)
       )
fig.savefig('../plots/bulk_ReHeat')

# Contrasts
for contrast in np.unique(bf.contrast):
    x_label, y_label, s_label, c_label = get_contrast_corr(df, bf[bf.contrast==contrast])
    fig = dotplot(title='scHF vs Bulk ({0})'.format(contrast), 
            x=x_label, 
            y=y_label, 
            c=c_label, 
            s=s_label, 
            size_title='Overlap genes', 
            color_title='Spearman corr',
            num=0.75,
            figsize=(9,4)
           )
    fig.savefig('../plots/bulk_{0}'.format(contrast))