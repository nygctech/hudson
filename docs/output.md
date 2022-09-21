---
layout: default
title: Output
nav_order: 6
---


# Output

<p align="justify ">

hudson currently creates 2 outputs, with both of them being dictionaries. The output is pickled the output directory specified by the user in the pipleline config file. 
</p> 

## Cell Type Composition Dict

<p align="justify ">

This is the main dictionary, which has all the neighborhood compositions for each cell identified in the tissue. The keys of the dictionary are the node numbers and the corresponding values are the percentange cell type compositions of the their 2nd order neighborhood, stored in a pandas dataframe format.
</p> 


## Net Neighborhood Dict 

<p align="justify ">

This is the net neighborhood dictionary, with nodes and their 2nd order neighbors stored as they dictionary keys and values. 
</p> 
