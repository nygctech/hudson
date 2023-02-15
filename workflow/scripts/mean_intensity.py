import zarr
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client
import torch
import joblib
#from dask_jobqueue import SLURMCluster
import skimage
import time
import pickle
from os.path import exists, join
from joblib import Parallel, delayed
from joblib import parallel_backend
from dask.distributed import progress
from pathlib import Path
from utils import open_zarr, get_cluster, get_logger
from skimage.measure import regionprops_table



# open xarray image from zarr store
image = open_zarr(snakemake.input[0])

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])
logger.info(f'Opened {image.name} zarr')

# Open instance labels
labels = skimage.io.imread(Path(snakemake.input[1]))
logger.info(f'Opened {image.name} labels')


# Get markers
marker_list = list(image.marker.values)
for m in snakemake.config.get('mean_intensity',{}).get('exclude', []):
    marker_list.remove(m)
logger.info('Features:')
for m in marker_list:
    logger.info(m)

# plane_dict = {}

# for mark in marker_list:
#     try:
#         plane_dict.update({mark: image.sel(marker = mark)})
#     except:
#         logger.info(f'Could not find {mark} in image')
#         pass
    
if torch.cuda.is_available() == False:
    
#     label_props = ('label', 'area', 'centroid', 'eccentricity', 'equivalent_diameter_area', 'feret_diameter_max', 'orientation',
#                    'perimeter', 'perimeter_crofton')
#     cell_props = regionprops_table(labels, properties = label_props)
    
    mean_intensity_per_marker = {}
    for m in image.marker.values:
        logger.info(f'Measuring {m}')
        props = regionprops_table(labels, intensity_image = image.sel(marker = m).values, 
                                  properties = ('intensity_mean',))
        mean_intensity_per_marker.update({m:props['intensity_mean']})





    # Start dask cluster
    # specify default worker options in ~/.config/dask/jobqueue.yaml
#     winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
#     cluster = get_cluster(**winfo)
#     logger.debug(cluster.new_worker_spec())
#     logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
#     nworkers = snakemake.params.tiles*2
#     logger.info(f'Scale dask cluster to {nworkers}')
#     cluster.scale(nworkers)
#     client = Client(cluster)

#     val = np.max(labels)
    
    
#     def get_pixels(lab,pl):
#         m = plane_dict[pl].values[labels == lab+1].mean()
#         return m
    
#     mean_intensity_per_marker = {}
#     for plane in plane_dict.keys():
    
#         with parallel_backend('dask',scheduler_host=cluster.scheduler._address,wait_for_workers_timeout=20):
#             mean_int = Parallel(n_jobs=-1)(delayed(get_pixels)(lab, pl = plane) for lab in range(val))
#         mean_intensity_per_marker.update({plane:mean_int})
    
    logger.info(f'Writing features')
    with open(Path(snakemake.output[0]), 'wb') as f:
        pickle.dump(mean_intensity_per_marker, f)
        
#    client.close()
#    cluster.close()
        
    
else:
    
    plane_dict = {}
    for mark in marker_list:
        try:
            plane_dict.update({mark: image.sel(marker = mark)})
        except:
            logger.info(f'Could not find {mark} in image')
            pass

    def get_mean_intensity(pl):
        result_ar = np.zeros(mx)
        tr = torch.from_numpy(pl)
        for r in range(mx):
            result_ar[r] = (tr[lab == r+1]).float().mean()
        return result_ar

    lab = torch.from_numpy(labels.astype('int'))
    mx = np.max(labels)
    
    mean_intensity_per_marker = {}
    for plane in plane_dict.keys():
        pl = plane_dict[plane].values
        mean_int = get_mean_intensity(pl)
        mean_intensity_per_marker.update({plane:mean_int})
    
    with open(Path(snakemake.output[0]), 'wb') as f:
        pickle.dump(mean_intensity_per_marker, f)
    
        

    
