---
layout: default
title: Overview
nav_order: 1
---


## Overview


hudson is a computational pipeline for spatial analysis of single cells in data obtained via multiplexed fluorescence imaging. It is developed by *Jagjit Singh, Kunal Pandit, and Rui Fu*. The pipeline can installed by following the steps [here](https://github.com/nygctech/hudson)

Fluorescent Imaging, is the visualization of fluorescent dyes or proteins as labels for molecular processes or structures. It enables a wide range of experimental observations including the location and dynamics of gene expression, protein expression and molecular interactions in cells and tissues. fluorescence microscopy. 

Following is an example of a multiplexed fluorescence image:
![Image](https://user-images.githubusercontent.com/42875353/185256327-27dfeb89-2cce-4bb7-b617-a434e7cf65dd.png){:height="50%" width="50%"}. 


The pipeline has parallel capibility meaning it can be scaled over any cluster using the SLURM manager. It makes use of GPU nodes if avaiable for quick image processing and computation. The structure is defined in the next section.  
