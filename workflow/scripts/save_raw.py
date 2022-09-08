from pre.utils import get_config
from pre import image_analysis as ia
from utils import get_cluster, get_logger
from dask.distributed import Client, wait
from pathlib import Path
import yaml
from math import ceil

# Find raw image path
experiment_config = get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])
image_path = Path(exp_dir) / Path(image_path)

# Get section name
section_name = snakemake.params.section

# Start logger
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])

# Open raw images
image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name,
                           logname = f'{section_name}.image')

# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
ntiles = int(len(image.im.col)/2048)
min_workers = max(1,2*ntiles)
max_workers = 4*min_workers

# Print out info about section
logger.info(f'machine:: {image.machine}')
logger.info(f'image path:: {image_path}')
logger.info(f'section:: {section_name}')

# Write raw images as zarr store
with Client(cluster) as client:

    cluster.adapt(minimum = min_workers, maximum=max_workers)
    client.wait_for_workers(min_workers, 60*10)


    # Write Raw Images
    save_path = Path(snakemake.output[0]).parents[0]
    delayed_store = image.save_zarr(save_path, compute = False)
    logger.info('Writing Raw Images')
    future_store = client.persist(delayed_store)
    wait(future_store)
    logger.info('Finished Writing Raw Images')
        
# write section info to file
section_info = {'tiles': ntiles,
                'planesizeMB':ceil(image.im.nbytes/(1024*1024)),
                'rawpath': str(save_path),
                'machine': image.machine,
                'experiment': experiment_config['experiment']['experiment name']
               }

with open(snakemake.output[1], 'w') as f:
    f.write(yaml.dump(section_info))
    