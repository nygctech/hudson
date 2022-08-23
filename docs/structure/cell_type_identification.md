---
layout: default
title: Cell Type Identification
parent: Structure 
nav_order: 5
---
## Clustering based on protein intensity and cell spatial parameters

<p align="justify ">
We treat the protein intensity measurements and scalar spatial parameters for each segmented cell as separate modes of data, and use Seurat Weighed Nearest Neighbor Analysis (WNN) to cluster cells using both.
<p align="justify ">
  
## Cell type classification from reference single cell expression dataset
  
<p align="justify ">
Cell type classification is automated, comparing all protein intensities with a reference dataset. This can be previously labeled scRNA-seq or PySeq2500 4i data, stored as Seurat object in RDS format. Seurat Canonical Correlation Analysis (CCA) has been shown to be performant in cross-modal reference mapping.
<p align="justify ">
