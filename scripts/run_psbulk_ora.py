import pandas as pd
import numpy as np

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

# Read DEG
input_path = '../plot_data/deg/deg.csv'
df = pd.read_csv(input_path)

# Get gene set db
db = get_db()

# Get different conditions and cell types
conditions = np.unique(df['condition'])

# Run ORA for each condition and cell type
dfs = []
for condition in conditions:
    cell_types = np.unique(df[df['condition']==condition]['cell_type'])
    for cell_type in cell_types:
        # Filter by cond, ctype
        msk = (df['condition'] == condition) & (df['cell_type'] == cell_type) & (df['pvals'] < 0.05)
        deg = df[msk]['names'].tolist()
        # Run ORA and filter by adj pval
        df_ora = ORA(deg, db).sort_values('adj_pvalues')
        df_ora['cell_type'] = cell_type
        df_ora['condition'] = condition
        df_ora = df_ora[df_ora['adj_pvalues'] < 0.05]
        dfs.append(df_ora)
        
# Merge dfs into one
df_ora = pd.concat(dfs, ignore_index=True)

# Save
df_ora.to_csv('../plot_data/deg/ora.csv', index=False)
