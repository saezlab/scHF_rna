# scHF_scanpy

This is a sc pipeline based on `scanpy` designed to analyse a collection of samples from Hearth Failure (HF) patients.


Get diss genes:
```
mkdir data
cd data
wget https://raw.githubusercontent.com/kieranrcampbell/scrnaseq-digestion-paper/master/data/deliverables/coregene_df-FALSE-v3.csv
```

Steps:

run qc with thresholds
run qc plots

run integration
run integration plots

run clustering
run clus plotting

Packages:
```
scanpy
scrublet
harmonypy
```
