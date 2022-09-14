---
layout: default
title: Unmix
parent: Structure
nav_order: 4
---

# Unmix
Spectral spillover in images are removed using the unmixing parameters computed with PICASSOnn. Images are saved as marker instead as channel/cycle.

<xarray.DataArray 'm387ntga1' (marker: 13, obj_step: 10, row: 10532, col: 8960)>
dask.array<getitem, shape=(13, 10, 10532, 8960), dtype=uint16, chunksize=(1, 1, 10532, 1792), chunktype=numpy.ndarray>
Coordinates:
    channel   (marker) int64 dask.array<chunksize=(4,), meta=np.ndarray>
  * col       (col) int64 0 1 2 3 4 5 6 7 ... 8953 8954 8955 8956 8957 8958 8959
    cycle     (marker) int64 dask.array<chunksize=(4,), meta=np.ndarray>
  * marker    (marker) <U7 'LMNB1_1' 'ELAVL2' 'GFAP' ... 'PDGFRA' 'MBP'
  * obj_step  (obj_step) int64 8847 9082 9317 9552 ... 10257 10492 10727 10962
  * row       (row) int64 0 1 2 3 4 5 6 ... 10526 10527 10528 10529 10530 10531 

## Input
Zarr store of processed images
YAML unmixing parameters

## Output
Zarr store of final images

# Logs
Script log
Dask performance report
