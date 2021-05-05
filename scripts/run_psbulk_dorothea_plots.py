import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from plotting import dotplot


# Read dorothea results
results = pd.read_csv('../plot_data/func/dorothea.csv').sort_values(['names'])

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
              edgecolor=['black' if pval < 0.05 else 'white' for pval in df['pvals']],
              figsize=(18,6)
             )
    
    # Save
    fig.savefig('../plots/psbulk_dorothea_{0}.png'.format(contrast), bbox_inches='tight')