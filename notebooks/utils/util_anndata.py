#!/usr/bin/env python
# coding: utf-8

import dask.dataframe as dd
import dask.array as da
import dask.bag as db
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import skimage
import anndata as ad
import pandas as pd


def return_anndata(path):
    
    img = skimage.io.imread(path)
    ski_img = skimage.measure.regionprops(img)
    df = pd.DataFrame(ski_img)
    var_array = df.head(1).values
    list_var = list(var_array.flatten())
    lables = len(ski_img)
    df = pd.DataFrame(index = list_var, columns = np.arange(1,lables,1))
    
    for ind in df.index:
        for val in df.columns:
            df.loc[ind,val] = getattr(ski_img[val],ind)
            
    for val in df.columns:
        df.loc['bbox',val] = np.asarray(df.loc['bbox',val]).reshape(1,4)
        
        
    for val in df.columns:
        df.loc['centroid',val] = np.asarray(df.loc['centroid',val]).reshape(1,2)
        
        
    for ind in df.index:
        for val in df.columns:
            if type(df.loc[ind,val]) == np.ndarray:
                df.loc[ind,val] = df.loc[ind,val][~np.isnan(df.loc[ind,val])].flatten()
                
                
                
    scalar_list = ['area','area_bbox','area_convex','area_filled','axis_major_length','axis_minor_length',
            'eccentricity', 'equivalent_diameter_area','euler_number','extent','feret_diameter_max',
              'label','orientation','perimeter','perimeter_crofton','solidity'] 
    scalar_df = df.loc[['area','area_bbox','area_convex','area_filled','axis_major_length','axis_minor_length',
            'eccentricity', 'equivalent_diameter_area','euler_number','extent','feret_diameter_max',
              'label','orientation','perimeter','perimeter_crofton','solidity']]
                   
    multidim_list = ['bbox','centroid','centroid_local', 'coords', 'image', 'image_convex', 'image_filled',
                'inertia_tensor', 'inertia_tensor_eigvals', 'moments', 'moments_central' ,'moments_hu',
                'moments_normalized']
    multidim_df = df.loc[['bbox','centroid','centroid_local', 'coords', 'image', 'image_convex', 'image_filled',
                'inertia_tensor', 'inertia_tensor_eigvals', 'moments', 'moments_central' ,'moments_hu',
                'moments_normalized']]
    
    
    mt = multidim_df.transpose()
    
    X = scalar_df.T.values
    adata = ad.AnnData(X)
    adata.var_names = scalar_df.index
    adata.obs_names = scalar_df.columns
    adata.var = pd.DataFrame(index = scalar_list, columns = scalar_list)
    adata.obs = pd.DataFrame(index = scalar_df.columns, columns = scalar_df.columns)
    for col in mt.columns:
        adata.obsm[col] = mt[col].values
        
    i = 0
    for row in multidim_list:
        adata.uns[row] = multidim_df.iloc[[i]].values.flatten()
        i = i+1
        
    return adata
    
        
        

    
    
    
    
                
                
    
        
    
            
    
                            
    
    

