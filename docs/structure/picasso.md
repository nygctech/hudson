---
layout: default
title: PICASSO
parent: Structure
nav_order: 3
---

# PICASSO

## Unmix spectral spillover

![](https://user-images.githubusercontent.com/72306584/176486552-50e1bca9-65fd-4466-8c92-a114e48d2278.gif)

## PICASSOnn

PICASSOnn is an algorithm to remove spillover fluorescence by minimizing the mutual information (MI) between sink and source images<sup>1</sup>. MI between pairs of sink and source images is estimated and minimized with a neural network (MINE<sup>2</sup>)using stochastic gradient descent and GPU acceleration.

## Mixing model

$$ unmixed sink = sink - \sum_{i} \alpha_i(source_i - \beta_i) $$

## References

1. Seo, J. et al. PICASSO allows ultra-multiplexed fluorescence imaging of spatially overlapping proteins without reference spectra measurements. Nat Commun 13, 2475 (2022).
2. Belghazi, M. I. et al. MINE: Mutual Information Neural Estimation. arXiv:1801.04062 [cs, stat] (2018).

## Input
Zarr store of processed images
YAML section summary

## Output
YAML unmixing parameters

## Log
Script log
