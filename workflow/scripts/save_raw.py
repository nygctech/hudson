from pre.utils import get_config
from pre import image_analysis as ia
from utils import get_cluster, get_logger, position
from dask.distributed import Client, wait, performance_report
from pathlib import Path
import yaml
from math import ceil
from shutil import copy2
from os import makedirs
from glob import iglob

# Find raw image path
experiment_config = get_config(snakemake.input[0])
exp_dir = Path(snakemake.config['experiment_directory'])
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])
image_path = exp_dir / image_path

# Get section name
section_name = snakemake.params.section

# Start logger
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])

# Open raw images
image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name,
                           logname = f'{section_name}.image')
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
logger.info(f'section:: {section_name}')



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
        
    # Copy Experiment Config     
    log_dir = Path(snakemake.output[2])
    os.makedirs(log_dir, exist_ok = True) 
    if not (log_dir / 'config.cfg').exists():
        copy2(snakemake.input[0], log_dir / 'config.cfg')

    # Copy Experiment and flowcell logs
    for f in iglob(str(exp_dir / 'logs' /'*.log'):
        if not (log_dir / f.name).exists():
            copy2(f, log_dir / f.name)

    # Copy over autofocus data for the section
    focus_config = exp_dir / 'logs' / exp_name / 'focus_config.cfg'
    if not (log_dir / 'focus_config.cfg').exists():
	copy2(focus_config, exp_dir / 'focus_config.cfg')
    focus_dir = Path(snakemake.output[3])
    os.makedirs(focus_dir, exist_ok = True)
    pos_dict = position(experiment_config.get('section', section_name))
    section_id = str(pos_dict['x_center'])+str(pos_dict['y_center'])
    for f in iglob(str(focus_dir / section_id + '*')):
        copy2(f, focus_dir / f.name)
