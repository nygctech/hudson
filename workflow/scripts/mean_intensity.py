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

plane = image.sel(marker = 'LMN1b')


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
                    #extra=["--lifetime", "55m", "--lifetime-stagger", "4m"])
        client = Client(cluster, timeout="50s")

        return cluster, client

    cluster, client = get_cluster()

    def scale_cluster(count): 
        cluster.scale(count)
        return cluster.dashboard_link
    scale_cluster(5)

    #Way 3: Computing the same using dask compute without persisting Data
    val = np.max(labels)
    def get_pixels(lab):
        m = plane.values[labels == lab+1].mean()
        return m

    with parallel_backend('dask',scheduler_host=cluster.scheduler._address):
        results = Parallel(n_jobs=-1)(delayed(get_pixels)(lab) for lab in range(val))
        
    print(results)
    
    intensity_dict = {}
    intensity_dict.update({'values':results})
    with open(Path(snakemake.output[0], 'wb')) as f:
        pickle.dump(intensity_dict, f)
        
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
    pl = plane.values
    mean_int = get_mean_intensity(pl)
    
    intensity_dict = {}
    intensity_dict.update({'values':mean_int})
    with open(snakemake.output[0], 'wb') as f:
        pickle.dump(intensity_dict, f)
    
        

    