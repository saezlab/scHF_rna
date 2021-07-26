import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from plotting import varshift_condition, varshift_cell_type 
    
# Open data
df = pd.read_csv('../plot_data/var_shifts/shifts.csv')

# By condition
varshift_condition(df[df.type == 'rna'], 'psbulk_varshifts_rna_condition.png')
varshift_condition(df[df.type == 'atac'], 'psbulk_varshifts_atac_condition.png')
varshift_condition(df, 'psbulk_varshifts_both_condition.png')
    
# By cell type
varshift_cell_type(df[df.type == 'rna'], 'psbulk_varshifts_rna_celltype.png')
varshift_cell_type(df[df.type == 'atac'], 'psbulk_varshifts_atac_celltype.png')    
varshift_cell_type(df, 'psbulk_varshifts_both_celltype.png')