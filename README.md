# Bulk RNA-seq analysis on PDAC tissu : SCIE-PANC project (Vera Pancaldi Lab, CRCT)

This repository contains all the work and analyis done on bulk RNA-seq data. The primary objective of this project was to explore the data and to connect clinical data (for example : response to treatment) and the tumoral micro-environment. The analysis includes PCA, sPLS-DA, heatmap, deconvolution with multideconv tool, and cell niches (CellTFusion tool)

## Dataset
248 samples from patients with PDAC, tumour tissue from biopsies of primary or metastatic tumours. Bulk RNA-seq data, STAR v2.7.4a alignment to the human genome hg38.

Clinical data for all the samples (mutation type, response to therapy, moffitt...) 

## Resources

R and RStudio (version ≥ 4.3.x)

Additional dependencies listed in the renv environment (see the env folder)

The project uses renv to manage R dependencies, ensuring reproducibility. The environment is compatible with R version ≥ 4.3.x

## Author
Augustine Blanc-Boekholt

Vera Pancaldi Lab, NetB(IO)2, CRCT 
