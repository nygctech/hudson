---
layout: default
title: Micro Environment 
parent: Structure
nav_order: 10
---

## Micro Environment Compositions

<p align="justify ">

Now, that all our essential data is condensed into a single object, we extract the centroids from our anndata object and using these centroids, we create a voronoi tessalation/diagram. More information on Voronoi transformations can be found here. The purpose of the voronoi tessellation is to form localized regions over tissue section. This lets us build micro-environments for localized spatial analysis of cell type compositions. Hudson computes the cell type composition within each micro-environment and stores it for downstream analysis. 
  
</p>

