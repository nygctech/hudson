---
layout: default
title: Overview
nav_order: 1
---


## Overview

            

<p align="justify ">
 Hudson is a computational pipeline for spatial analysis of tissue sections from multiplexed fluorescence imaging obtained from a converted Illumina
 HiSeq2500 sequencing system. It is developed by <b> Jagjit Singh, Kunal Pandit, Rui Fu, and Sanja Vickovic. </b>
</p>
 
<p align="justify ">
 Please follow the steps on the left for an in dept view and understanding of hudson. This documentationalso includes examples for the user to follow and test on their own system. The preprint is available here
</p>

<p align="justify ">
 Following is an example of a multiplexed fluorescence image:
</p>


<p align="center">
  <img width="50%" height="50%" src="https://github.com/nygctech/hudson/blob/docs/docs/spinal_tissue.png">
</p>


<p align="justify ">
 Analysing such images is computationally intensive, thus the pipeline has parallel compute capability. It can be scaled over any computer cluster using
 the Slurm manager. Slurm is an open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux
 clusters. Additionally, It makes use of GPUs or GPU nodes if avaiable for even quicker image processing and computation. The structure is defined in the
 next section.
</p>

  ![Image](banner.png)


