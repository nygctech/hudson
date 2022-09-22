---
layout: default
title: Spatial Neighborhood
parent: Structure
nav_order: 9
---

# Spatial Neighborhood

## Voronoi Diagram

<p align="justify ">

Now, that all our essential data is condensed into a single object, we extract the centroids from our anndata object and using these centroids, we create a voronoi tessalation/diagram. More information on Voronoi transformations can be found here. The purpose of the voronoi tessellation is to form localized regions over tissue section. This lets us build micro-environments for localized spatial analysis of cell type compositions. Hudson computes the cell type composition within each micro-environment and stores it for downstream analysis. 
</p>


## Graph Network

<p align="justify ">

Hudson also builds a network of cell connectivities using a K-Nearest-Neighbor approach. This graph can then be used for local spatial autocorrelation analysis based nuclei or cell based markers. The pipline stores the spatial autocorrelation coefficent (Moran's I) for each marker available. 
It also stores the KNN graph network and the user can fetch neighboring connections for any cell and neighborhood order. Furthermore, the 2nd order neighborhood connections are pickled as a dictionary, for quick analysis. 
 
</p>
