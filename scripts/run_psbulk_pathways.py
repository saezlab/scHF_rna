import numpy as np
import pandas as pd

import os

from utils import run_mlm

df = pd.read_csv('../plot_data/deg/psbulk_deg.csv')
df['index'] = df['contrast'] + '$' + df['cell_type']
contrasts = ['HF-AvsHealthy', 'HF-CKDvsHealthy','HF-CKDvsHF-A']

mat = df.pivot(index='names', columns='index', values='tvals')
mat[np.isnan(mat)] = 0

network = pd.read_csv('../data/progeny_100.csv')

# Run mlm
acts = run_mlm(mat, network)
acts[['contrast', 'cell_type']] = acts['condition'].str.split('$', 1, expand=True)
acts = acts[['source', 'contrast', 'cell_type', 'score', 'p_value']]
acts = acts.rename(columns={'score' : 'tvals'})

os.makedirs('../plot_data/func/', exist_ok=True)
acts.to_csv('../plot_data/func/pathways.csv', index=False)