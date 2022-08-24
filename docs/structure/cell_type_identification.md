---
layout: default
title: Cell Type Identification
parent: Structure 
nav_order: 5
---
# Cell Type Identification

## Preprocessing and clustering based on protein intensity and cell spatial parameters

<p align="justify">
We treat the protein intensity measurements and scalar spatial parameters for each segmented cell as separate modes of data to be separatly normalized, scaled, and batch-corrected. Then, Seurat Weighed Nearest Neighbor Analysis (WNN) is used to cluster cells and project into two-dimensional space using both modalities. This also allows for flexibility in integrating multiple datasets.
<p/>
  
## Cell type classification from reference single cell expression dataset
  
<p align="justify">
Automated cell type classification is achieved through comparing all protein intensities with a reference dataset. This can be previously labeled scRNA-seq or PySeq2500 4i data, stored as Seurat object in RDS format, as specified in the experiment `config/config.yaml` file (and metadata column name that contains cell type information in the reference). The underlying algorithm used here, Seurat Canonical Correlation Analysis (CCA), has been shown to be performant in cross-modal reference mapping. Cell type results are saved as a csv file in the output `tables` directory.
<p/>
