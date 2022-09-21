---
layout: default
title: Intensity Computation
parent: Structure
nav_order: 6
---
# Intensity Computation 

<p align="justify ">
  This step in the pipeline computes mean intensities for all available markers in each label. This is one the most computationally intensive steps is
  crucial for detecting cell type composition. Following is an example of the mean pixel intensity per label in the case of Lamin-B1 and MBP nucelic marker.
</p> 

`{'LMN1b': [374.4911504424779,564.6321839080459,...] 'MBP': [134.27433628318585,181.33333333333334,..]}`

<p align="justify ">
  As we can see ablove the pixel intensity per marker is saved as a 1d array for all available markers in the image. All the intensity information saved as
  a dictionary with the markers being keys and the values being the corresponding arrays containg the mean intensities. 
</p> 

<p align="justify ">
  The pipeline will automatically scale the code over the compute cluster. The computation is further accelarated if a Graphical Processing Unit is present
  in the system or cluster. 
</p> 
