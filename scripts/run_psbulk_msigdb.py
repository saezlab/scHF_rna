import numpy as np
import pandas as pd

import os

from utils import run_ora, get_db

df = pd.read_csv('../plot_data/deg/psbulk_deg.csv')
df['index'] = df['contrast'] + '$' + df['cell_type']
contrasts = ['HF-AvsHealthy', 'HF-CKDvsHealthy','HF-CKDvsHF-A']

mat = df.pivot(index='names', columns='index', values='tvals')
mat[np.isnan(mat)] = 0

network = []
for row in get_db():
    gset_id, gset_desc, gset = row
    for gene in gset:
        network.append([gset_id, gene, 1, 1])
network = pd.DataFrame(network, columns=['source', 'target', 'mor', 'likelihood'])
network = network.drop_duplicates(['source', 'target'])

# Run ora
acts = run_ora(mat, network)
acts[['contrast', 'cell_type']] = acts['condition'].str.split('$', 1, expand=True)
acts = acts[['source', 'contrast', 'cell_type', 'score', 'p_value']]
acts = acts.rename(columns={'score' : 'tvals'})

os.makedirs('../plot_data/func/', exist_ok=True)
acts.to_csv('../plot_data/func/msigdb.csv', index=False)
