import zarr
import xarray as xr
import numpy as np
import numpy as np
import dask
from pyseq import image_analysis as ia
from dask.distributed import Client
import torch
import joblib
from dask_jobqueue import SLURMCluster
import skimage
import time
import pickle
from os.path import exists, join
from joblib import Parallel, delayed
from joblib import parallel_backend
from dask.distributed import progress
from pathlib import Path


labels = skimage.io.imread(Path(snakemake.input[1]))
image_path = Path(snakemake.input[0])
im_name = image_path.stem
image = xr.open_zarr(image_path).to_array()
image = image.squeeze().drop_vars('variable').rename(im_name)

marker_list = ['LMN1b', 'GFAP','ELAVL2','MBP','PVALB']

plane_dict = {}

for mark in marker_list:
    try:
        plane_dict.update({mark: image.sel(marker = mark)})
    except:
        pass
    
if torch.cuda.is_available() == False:

    def get_cluster(queue_name = 'pe2', log_dir=None):
        """ Make dask cluster w/ workers = 2 cores, 32 G mem, and 1 hr wall time.

            return cluster, client
        """

        cluster = SLURMCluster(
                    queue = queue_name, 
                    cores = 6 ,
                    memory = '60G',
                    walltime='1:00:00')
                    
        client = Client(cluster, timeout="50s")

        return cluster, client

    cluster, client = get_cluster()

    def scale_cluster(count): 
        cluster.scale(count)
        return cluster.dashboard_link
    scale_cluster(5)

    
    val = np.max(labels)
    
    
    def get_pixels(lab,pl):
        m = plane_dict[pl].values[labels == lab+1].mean()
        return m
    
    mean_intensity_per_marker = {}
    for plane in plane_dict.keys():
    
        with parallel_backend('dask',scheduler_host=cluster.scheduler._address,wait_for_workers_timeout=20):
            mean_int = Parallel(n_jobs=-1)(delayed(get_pixels)(lab, pl = plane) for lab in range(val))
        mean_intensity_per_marker.update({plane:mean_int})
    

    with open(Path(snakemake.output[0], 'wb')) as f:
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
    
    with open(Path(snakemake.output[0], 'wb')) as f:
        pickle.dump(mean_intensity_per_marker, f)
    
        

    