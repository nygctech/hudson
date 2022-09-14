---
layout: default
title: Segmentation
parent: Structure
nav_order: 5
---

# Segmentation

<p align="justify ">

 In the semgmenation step, the pipline performs instance segmentation on the image to generate a masks array of the same dimension as the input. The masks array is essentially the assignment of each pixel to a certain cell in the image. The users can select over which marker(LMN1b, GFAP, MBP, EVALV2,..) they want to segment on. Each cell is indexed as a label starting from 1 to the n'th cell that is segmented.

</p>

Below, on the left is the original images and on the right is the generated masks.
