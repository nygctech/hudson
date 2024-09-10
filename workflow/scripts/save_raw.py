from pre.utils import get_config
from pre import image_analysis as ia
from utils import get_cluster, get_logger
from dask.distributed import Client, wait, performance_report
from pathlib import Path
import yaml
from math import ceil

# Find raw image path
experiment_config = get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])


# Get section name
section_name = snakemake.params.section

# image_path = Path(image_path) / Path(section_name)
image_path = Path(image_path)
# Start logger
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])
logger.info(f'path:: {image_path}')
logger.info(f'section:: {section_name}')

# Open raw images
image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name,
                           logname = f'{section_name}.image')
logger.info(f'machine::{image.machine}')
logger.debug(f'{image.im.shape}')

assert image.machine != ''

# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
ntiles = image.im.col.size//2048
nworkers = max(1,ntiles)
logger.info(f'Scale dask cluster to {nworkers}')
cluster.scale(nworkers)
client = Client(cluster)

# Print out info about section
logger.info(f'machine:: {image.machine}')
logger.info(f'image path:: {image_path}')




# Write Raw Images
save_path = Path(snakemake.output[0]).parents[0]
delayed_store = image.save_zarr(save_path, compute = False)


# Send write to cluster
logger.info('Writing Raw Images')
with performance_report(filename=snakemake.log[1]):
    future_store = client.persist(delayed_store, retries = 10)
    futures = list(future_store.dask.values())
    wait(future_store)

    
# Double check no errors
futures_done = [f.done() for f in futures]    
if all(futures_done):
    logger.info('Finished Writing Raw Images')
        
    # get 1 plane and number of cycles
    sel = {}
    cycles = 1
    for key, value in  image.im.coords.items():
        sel[key] = value[0]
        if 'cycle' == key:
            cycles = len(value)
    plane = image.im.sel(sel) 
    
    
        
    # write section info to file
    section_info = {'tiles': ntiles,
                    'planesizeMB':ceil(plane.nbytes/(1024*1024)),
                    'rawpath': str(save_path),
                    'machine': image.machine,
                    'experiment': experiment_config['experiment']['experiment name'],
                    'cycles': cycles
                   }
    with open(snakemake.output[1], 'w') as f:
        f.write(yaml.dump(section_info))
        
   
    with open(snakemake.output[2], 'w') as f:
        experiment_config.write(f)
