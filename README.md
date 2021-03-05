# scHF_scanpy

This is a sc pipeline based on `scanpy` designed to analyse a collection of samples from Hearth Failure (HF) patients.

## Set up
Required packages:
```
scanpy
scrublet
bbknn
leidenalg
scCODA
seaborn
gseapy
```

Get dissociation genes:
```
mkdir data
cd data
wget https://raw.githubusercontent.com/kieranrcampbell/scrnaseq-digestion-paper/master/data/deliverables/coregene_df-FALSE-v3.csv
```
## Process data
### QC
To run qc filtering simply type, specific thresholds for `scrublet` can be set up per sample:
```
cd scripts/
python run_qc.py
```

Then to generate plots run:
```
python run_qc_plots.py
```
A QC plot will be generated per sample, example:
![Example QC single sample](/plots/qc_CK128.png)
And a summary of all samples:
![Summary QC](/plots/qc_summary.png)

### Sample integration
Steps:
1. Find HVG per batch
2. Filter top 3000 HVG that are HV in the maximum number of samples
3. Run PCA
4. Run `BBKNN`
5. Run UMAP

To run the integration, type:
```
python run_integration.py
```
Then this to plot the results
```
python run_integration_plots.py
```
These are the HVG QC plots:
![Summary QC HVG](/plots/hvg.png)

These are the integration projections:
![Integration projections](/plots/proj_integration.png)

### Clustering
We run a sequence of different leiden resolutions by:
```
python run_clus.py
```
Then
```
python run_clus_plots.py
```

We get a projection plot for each resolution with a dotplot using some given markers:
![Example leiden resolution](/plots/leiden_res_1.0.png)

### Annotation
Once we check the different resolutions, we manually annotate the clusters by running:
```
python run_annot.py
python run_annot_plots.py
```
This is the final annotation:
![Annotation](/plots/proj_annotation.png)


