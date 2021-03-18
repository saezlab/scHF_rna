import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from plotting import plot_ora
    
    
# Read ORA df
input_path = '../plot_data/deg/ora.csv'
df = pd.read_csv(input_path)

# Get different conditions and cell types
contrasts = np.unique(df['contrast'])

for contrast in contrasts:
    ora = df[df['contrast'] == contrast]
    cell_types = np.unique(ora['cell_type'])
    fig, axes = plt.subplots(len(cell_types),1, dpi=150, figsize=(4, 4*len(cell_types)), sharex=True)
    fig.suptitle('ORA DEG {0}'.format(contrast), fontsize=11)
    axes = axes.flatten()
    for cell_type, ax in zip(cell_types, axes):
        plot_ora(cell_type, ora[ora['cell_type'] == cell_type], ax=ax)
    
    fig.subplots_adjust(top=0.95)
    fig.tight_layout()
    fig.set_facecolor('white')
    
    fig.savefig('../plots/psbulk_ora_{0}.png'.format(contrast), bbox_inches='tight')