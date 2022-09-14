---
layout: default
title: Fix Lighting
parent: Structure
nav_order: 2
---

# Fix Lighting

## Correct Background
Hamamatsu line scanning CCD Time Delay Integration (TDI) cameras are used in the Illumina HiSeq2500 sequencing system. The pixels in these cameras are arranged in a line in groups of 256 or 512 pixels (depending on the model). Each group has a different dark pixel value which is corrected to create a flat background.

![](https://user-images.githubusercontent.com/72306584/190162480-64e540c5-f015-4fd0-a98b-66af51e40e5a.png)

## Register Channels

The HiSeq images 4 emission channels simultaneously, however the images from these 4 channels are not perfectly aligned and differ from machine to machine.
Machine specific sub pixel channel shifts are precomputed during machine commissioning using phase cross correlation. Channel shifts are applied using a rigid affine transformation.

![](https://user-images.githubusercontent.com/72306584/190162640-5f1679e8-2dd7-4322-91e2-788f61364ae1.png)

## Stitch Tiles
Since the HiSeq uses line scanning TDI camera, image tiles are arbitrarily long in the row dimension (Y axis) and have a fixed column width (X axis) of 2048 pixel. Tiles can also optionally overlap a specified pixel width across the X axis during image acquistion. Tiles are stitched across the X axis with overlapping regions removed.

## Processed Image

![](https://user-images.githubusercontent.com/72306584/190162778-d62e346a-3308-4cc7-861c-08428ced4a2c.png)

## Input
Zarr store of raw images

## Output
Zarr store of processed images

## Logs
Script log
Dask performance report
