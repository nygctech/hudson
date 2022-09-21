---
layout: default
title: Data Object
parent: Structure
nav_order: 7
---

# Data Object

<p align="justify ">
  Highly multiplexed images and can have a lot of characteristic data. These are essentially the region proerties of the image which characterize it. More
  importantly, the number of labels.
</p>

<p align="justify ">
  
  In the step , the image labels along with their cell type composition and location is condensed and stored in an Anndata object. Anndata is a Python
  package for handling annotated data matrices in memory and on disk, positioned between pandas and xarray. All data is condensed to a 2 dimensional
  matrix, with columns being the labels and rows being the markers, the coordinates, and other relevant region properties. Below is the summary of the
  Anndata Object built from the example dataset. 

</p>

<p align="justify ">
  The regions properties are then stacked in a 2 dimensional matrix, along with the mean intensities per label for each of the protein markers, we covered
  in the mean intensity section. There is additional metadata and variable information added to object. Following condensed used gzip into an H5ad object.
  H5ad is anndata's native storage format which is based on hdf5.
</p>
