import pandas as pd
import numpy as np

from utils import get_db, ORA

# Read DEG
input_path = '../plot_data/deg/psbulk_deg.csv'
df = pd.read_csv(input_path)

# Get gene set db
db = get_db()

# Get different conditions and cell types
contrasts = np.unique(df['contrast'])

# Run ORA for each condition and cell type
dfs = []
for contrast in contrasts:
    cell_types = np.unique(df[df['contrast']==contrast]['cell_type'])
    for cell_type in cell_types:
        # Filter by cond, ctype
        msk = (df['contrast'] == contrast) & (df['cell_type'] == cell_type) & (df['pvals'] < 0.05)
        deg = df[msk]['names'].tolist()
        # Run ORA and filter by adj pval
        df_ora = ORA(deg, db).sort_values('adj_pvalues')
        df_ora['cell_type'] = cell_type
        df_ora['contrast'] = contrast
        df_ora = df_ora[df_ora['adj_pvalues'] < 0.05]
        dfs.append(df_ora)
        
# Merge dfs into one
df_ora = pd.concat(dfs, ignore_index=True)

# Save
df_ora.to_csv('../plot_data/deg/ora.csv', index=False)
