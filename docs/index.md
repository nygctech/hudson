---
layout: default
title: Overview
nav_order: 1
---


## Overview



 hudson is a computational pipeline for spatial analysis of single cells in data obtained via multiplexed fluorescence imaging. It is developed by **Jagjit Singh, Kunal Pandit, and Rui Fu**.
  
<p align="justify ">
 Fluorescent Imaging is the visualization of fluorescent dyes or proteins as labels for molecular processes or structures.It enables a wide range of 
 experimental observations including the location and dynamics of gene expression, protein expression and molecular interactions in cells and tissues.
</p> 
  
<p align="justify ">
 Multiplexed imaging is an emerging and exciting way to extract information from human tissue samples by visualizing many more biomarkers than traditional
 microscopy. By observing many biomarkers simultaneously, biological pathways previously explored only in isolation can be explored in concert, and
 complex tissue and cell phenotypes can be identified and probed. 
</p> 
  
 Following is an example of a multiplexed fluorescence image:
  

  ![Image](spinal_tissue.png){:height="50%" width="50%"}.
 
<p align="justify ">
  
  We can see that the image has different colours which correspond to different..
</p> 

<p align="justify ">
 
  Analysing such images is computationally intensive, thus the pipeline has parallel compute capibility. It can be scaled over any computer cluster using
  the Slurm manager. Slurm is an open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux
  clusters.
</p> 
  
  ![async-embarrassing](async_embarrassing_main.gif){:height="70%" width="70%"}.

<p align="justify ">
  Additionally, It makes use of GPUs or GPU nodes if avaiable for even quicker image processing and computation. The structure is defined in the next
  section. 
</p> 

