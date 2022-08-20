---
layout: default
title: Overview
nav_order: 1
---


## Overview


hudson is a computational pipeline for spatial analysis of single cells in data obtained via multiplexed fluorescence imaging. It is developed by **Jagjit Singh, Kunal Pandit, and Rui Fu**. 

Fluorescent Imaging is the visualization of fluorescent dyes or proteins as labels for molecular processes or structures. It enables a wide range of experimental observations including the location and dynamics of gene expression, protein expression and molecular interactions in cells and tissues.

Multiplexed imaging is an emerging and exciting way to extract information from human tissue samples by visualizing many more biomarkers than traditional microscopy. By observing many biomarkers simultaneously, biological pathways previously explored only in isolation can be explored in concert, and complex tissue and cell phenotypes can be identified and probed. 

Following is an example of a multiplexed fluorescence image:
![Image](https://user-images.githubusercontent.com/42875353/185256327-27dfeb89-2cce-4bb7-b617-a434e7cf65dd.png){:height="50%" width="50%"}. 


Analysing such images is computationally intensive thus pipeline has parallel compute capibility. It can be scaled over any cluster using the SLURM manager.Additionally, It makes use of GPUs or GPU nodes if avaiable for even quicker image processing and computation. The structure is defined in the next section.  
