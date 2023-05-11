import pickle 
from pathlib import Path

import numpy as np
import skimage
import anndata as ad
import pandas as pd
import squidpy as sq
import muon as mu



label_im_input = snakemake.input[0]
intensity_pkl = snakemake.input[1]

# Get Morphological Features
features = ('area','area_bbox','area_convex','area_filled','axis_major_length','axis_minor_length',
            'eccentricity', 'equivalent_diameter_area','euler_number','extent','feret_diameter_max',
            'orientation','perimeter','perimeter_crofton','solidity', 'label', 'centroid')
label_im = skimage.io.imread(label_im_input)
feature_table = skimage.measure.regionprops_table(label_im, properties=features)

# Organize Morphological Features
coords = np.vstack([feature_table['centroid-0'], feature_table['centroid-1']]).T
ind = feature_table['label']
del feature_table['centroid-0']
del feature_table['centroid-1']
del feature_table['label']

# Make Morphological AnnData object
morph_df = pd.DataFrame(data = feature_table, index = ind)
morph_ad = ad.AnnData(X = morph_df, obsm = {'spatial':coords})


# Add graph
radius = None # get from config file
if radius is None:
    radius = (feature_table['feret_diameter_max'].mean()+feature_table['feret_diameter_max'].std())*3/2
print(f'used radius {radius} px')
sq.gr.spatial_neighbors(morph_ad, n_neighs = True, delaunay=True, radius=(0,radius), coord_type = 'generic')


# Make Protein Intensities AnnData Object
with open(Path(intensity_pkl), 'rb') as f:
    protein_dict = pickle.load(f)
protein_df = pd.DataFrame(data = protein_dict, index = ind)
protein_ad = ad.AnnData(X = protein_df)     
          

# Combine into muon object
mdata = mu.MuData({'protein': protein_ad, 'morphological': morph_ad})
# h5 can't save radius parameter as tuple, so save as float
mdata['morphological'].uns['spatial_neighbors']['params']['radius'] = radius
mdata.write(snakemake.output[0])