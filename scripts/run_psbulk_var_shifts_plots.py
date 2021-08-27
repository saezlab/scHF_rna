import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from plotting import varshift_condition, varshift_cell_type, cond_colors
    
# Open data
df = pd.read_csv('../plot_data/var_shifts/shifts.csv')
tab_celltypes = pd.read_csv('../tables/shifts_celltypes.csv')
tab_celltypes = tab_celltypes.set_index(['condition','cell_type'])
tab_conditions = pd.read_csv('../tables/shifts_conditions.csv')

# By condition
#varshift_condition(df[df.type == 'rna'], 'psbulk_varshifts_rna_condition.svg')
#varshift_condition(df[df.type == 'atac'], 'psbulk_varshifts_atac_condition.svg')
#varshift_condition(df, 'psbulk_varshifts_both_condition.svg')
    
# By cell type
order = np.sort(df['cell_type'].unique())
hue_order = list(cond_colors.keys())
#varshift_cell_type(df[df.type == 'rna'], tab_celltypes, 'psbulk_varshifts_rna_celltype.svg')
#varshift_cell_type(df[df.type == 'atac'], tab_celltypes, 'psbulk_varshifts_atac_celltype.svg')    
varshift_cell_type(df, tab_celltypes, 'psbulk_varshifts_both_celltype.pdf', order, hue_order)
