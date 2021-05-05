import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from plotting import dotplot

# Read progeny results
results = pd.read_csv('../plot_data/func/progeny.csv').sort_values(['names', 'cell_type'])

# Get unique contrasts
contrasts = np.unique(results.contrast)

for contrast in contrasts:
    # Subset by contrast
    df = results[results.contrast == contrast]
    
    # Plot dot plot
    fig = dotplot(title=contrast,
              x=df['names'],
              y=df['cell_type'],
              c=df['logfoldchanges'],
              s=-np.log(df['pvals']),
              size_title='-log(pvalue)', 
              color_title='Mean change',
              edgecolor=['black' if pval < 0.05 else 'white' for pval in df['pvals']]
             )
    
    # Save
    fig.savefig('../plots/psbulk_progeny_{0}.png'.format(contrast), bbox_inches='tight')
