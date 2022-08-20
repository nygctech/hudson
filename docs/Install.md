---
layout: default
title: Installation
nav_order: 2
---

## Installing Hudson




The pipiline uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system to unite all its differnet components into one single workflow that can be replicated over any machine or scaled over any cluster. The different sections being: 

1) **Fix Lighting** : Explanation. <br>
2) **Unmix**: Explanation. <br>
3) **Segmentation** : In the semgmenation step, the pipline performs instance segmentation on the image to generate a masks array of the same dimension as the input. The masks array is essentially the assignment of each pixel to a certain cell in the image. The users can select over which marker(LMN1b, GFAP, MBP, EVALV2,..) they want to segment on. Each cell is indexed as a label starting from 1 to the n'th cell that is segmented. <br>
4) **Intensity Computation** : The next step is to compute the mean intensity per label using the masks array and the image. We need this to infer the cell types for the different labels that we got from the segmentation step. <br>
5) **Cell Type Indentification** <br>
6) **Annotation**: In the this we extract and condense the relevant image properties and cell type information in a single [anndata](https://anndata.readthedocs.io/en/latest/) object. This is very important to increase the efficieny, utility and scalability of the pipeline. <br>
7) **Spatial Neighborhood**:  Following this the pipeline will fetch the cell centroids from the Anndata object and perform a Voronoi tessellation. It will then use the Voronoi tessellation to create a graph of all cells using the edges and nodes from the tessellation. This graph will then used to build spatial neighborhoods, perform spatial analysis and store neighborhood connections for each cell. <br>
