---
layout: default
title: Data Object
parent: Structure
nav_order: 6
---

# Data Object

<p align="justify ">
Highly multiplexed images and can have a lot of characteristic data. These are essentially the region proerties of the image which characterize it. More importantly, the number of labels.
</p>

<p align="justify ">

Thus, it is important to store all the data in a single condensed memory efficient object. The pipeline makes use of Anndata, a data o.....library to meet the aforementioned needs. In this step hudson first computes the regions properties of the image, namely: _area, area_bbox, area_convex, area_filled, axis_major_length, axis_minor_length,eccentricity, equivalent_diameter_area euler_number, extent, feret_diameter_max, label, orientation, perimeter, perimeter_crofton','solidity'_ These properties are computed via the scikit-image libraries and in-detail definition can be found here.
</p>

<p align="justify ">

The regions properties are then stacked in a 2 dimensional matrix, along with the mean intensities per label for each of the protein markers, we covered in the mean intensity section. There is additional metadata and variable information added to object. Following condensed used gzip into an H5ad object. H5ad is anndata's native storage format which is based on hdf5.
</p>
