#!/usr/bin/env python
# coding: utf-8




import zarr
import numpy as np
import anndata
import scipy as sp
import dask 
from dask.distributed import Client
import xarray as xr
import pyseq
import skimage
from pyseq import image_analysis as ia
import torch
import pickle




labels = skimage.io.imread('/gpfs/commons/groups/nygcfaculty/PySeq/20210428_mouse_genotype_2/segmented_sections/m387ntga2_labels.tiff')
im = ia.get_HiSeqImages(image_path = 'zarrs/m387ntga2.zarr')

print(torch.cuda.is_available())


#format = one_z_plane_obj_step_channel_cycle
array_object_list = []
name_list = []
for i in im.im['channel'].values:
    for j in im.im['cycle'].values:
        #for k #in im.im['obj_step'].values[1]:
            k = 8028
            nme = "one_z_plane_"+str(i)+"_"+str(j)+"_"+str(k)
            name_list.append(nme)
            nme = im.im.sel(obj_step = k, cycle=j, channel = i)
            array_object_list.append(nme)







mx = labels.max()




def get_mean_intensity(pl):
    result_ar = np.zeros(mx)
    tr = torch.from_numpy(pl)
    for r in range(mx):
        result_ar[r] = (tr[lab == r+1]).float().mean()
    return result_ar




lab = torch.from_numpy(labels)
plane_mean_dict = {}
for i,nm in zip(range(len(array_object_list)),name_list):
    print(i,nm)
    pl = array_object_list[i].values
    mean_int = get_mean_intensity(pl)
    plane_mean_dict.update({nm:mean_int})
    with open('saved_dictionary_2.pkl', 'wb') as f:
        pickle.dump(plane_mean_dict, f)
        


    
    
    
    