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




# open xarray image from zarr store
image = open_zarr(snakemake.input[0])

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])
logger.info(f'Opened {image.name}')

# Open instance labels
labels = skimage.io.imread(Path(snakemake.input[1]))

# Get markers
marker_list = list(image.marker)
for m in snakemake.config['mean_intensity'].get('exclude', []):
    marker_list.remove(m)

plane_dict = {}

# Make sure only 1 objective step
if 'obj_step' in image.dims and 'obj_step' not in snakemake.config.get('segmentation',{}):
    if image.obj_step.size > 1:
        mid_step = image.obj_step[image.obj_step.size//2]
        image = image.sel(obj_step = mid_step)
        smk_logger.debug(image)


for mark in marker_list:
    try:
        plane_dict.update({mark: image.sel(marker = mark)})
    except:
        pass
    
if torch.cuda.is_available() == False:

    # Start dask cluster
    # specify default worker options in ~/.config/dask/jobqueue.yaml
    winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
    cluster = get_cluster(**winfo)
    logger.debug(cluster.new_worker_spec())
    logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
    nworkers = snakemake.params.tiles*2
    logger.info(f'Scale dask cluster to {nworkers}')
    cluster.scale(nworkers)
    client = Client(cluster)

    
    val = np.max(labels)
    
    
    def get_pixels(lab,pl):
        m = plane_dict[pl].values[labels == lab+1].mean()
        return m
    
    mean_intensity_per_marker = {}
    for plane in plane_dict.keys():
    
        with parallel_backend('dask',scheduler_host=cluster.scheduler._address,wait_for_workers_timeout=20):
            mean_int = Parallel(n_jobs=-1)(delayed(get_pixels)(lab, pl = plane) for lab in range(val))
        mean_intensity_per_marker.update({plane:mean_int})
    

    with open(Path(snakemake.output[0]), 'wb') as f:
        pickle.dump(mean_intensity_per_marker, f)
        
    client.close()
    cluster.close()
        
    
else:
    
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
    
        

    