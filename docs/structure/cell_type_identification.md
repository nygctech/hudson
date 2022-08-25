---
layout: default
title: Cell Type Identification
parent: Structure 
nav_order: 5
---
# Cell Type Identification
## Processing and clustering
<p align="justify">
We treat protein intensity measurements and scalar spatial parameters for each segmented cell as separate modes of data to be separatly normalized, scaled, and batch-corrected. Then, Seurat Weighed Nearest Neighbor Analysis <a href="https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html">(WNN)</a> is used to cluster cells and project into two-dimensional space using both modalities. This also allows for flexibility in integrating multiple datasets.
</p>
  
## Cell type classification from reference
<p align="justify">
Automated cell type classification is achieved through comparing all protein intensities with a reference dataset. This can be previously labeled scRNA-seq or PySeq2500 4i data, stored as Seurat object in RDS format, as specified in the experiment <code>config/config.yaml</code> file (and metadata column name that contains cell type information in the reference). The underlying algorithm used here, Seurat Canonical Correlation Analysis <a href="https://satijalab.org/seurat/reference/runcca">(CCA)</a>, has been shown to be performant in cross-modal reference mapping. Cell type results are saved as CSV file in the output <code>tables</code> directory, to be imported into anndata for downstream steps. If no appropriate reference is used, clustering results from above are saved instead.
</p>
