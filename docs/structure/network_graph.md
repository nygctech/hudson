---
layout: default
title: Netowork Graph
parent: Structure
nav_order: 9
---


## Cell Network Graph

<p align="justify ">

Hudson also builds a network of cell connectivities using a K-Nearest-Neighbor approach. This graph can then be used for local spatial autocorrelation analysis based nuclei or cell based markers. The pipline stores the spatial autocorrelation coefficent (Moran's I) for each marker available. 
It also stores the KNN graph network and the user can fetch neighboring connections for any cell and neighborhood order. Furthermore, the 2nd order neighborhood connections are pickled as a dictionary, for quick analysis. 
 
</p>
