---
layout: default
title: Save Raw
parent: Structure
nav_order: 1
---

# Save Raw
Save unprocessed images to a zarr store.

<xarray.DataArray 'm387ntga1' (channel: 4, cycle: 4, obj_step: 10, row: 10624, col: 12288)>
dask.array<getitem, shape=(4, 4, 10, 10624, 12288), dtype=uint16, chunksize=(1, 1, 1, 10624, 2048), chunktype=numpy.ndarray>
Coordinates:
  * channel   (channel) int64 558 610 687 740
  * cycle     (cycle) int64 1 2 3 4
  * obj_step  (obj_step) int64 8847 9082 9317 9552 ... 10257 10492 10727 10962
Dimensions without coordinates: row, col
Attributes:
    first_group:  0
    machine:      Origin
    scale:        1
    overlap:      256
    fixed_bg:     0

## Raw Image
![](https://user-images.githubusercontent.com/72306584/190162283-5f979f43-8b69-4920-ba01-2beddf13428d.png)

## Input
PySeq2500 experiment config

## Output
Zarr store of raw images
YAML section summary

## Logs
Script Log
Dask Performance Report
