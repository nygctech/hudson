from pre import image_analysis as ia
from utils import get_cluster, get_logger
from dask.distributed import Client, wait, performance_report
from pathlib import Path


# Get section name
section_name = Path(snakemake.input[0]).stem

# Start logger
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])

# Open image
image = ia.get_HiSeqImages(image_path = snakemake.input[0], logname = f'{section_name}.image')

# Check Raw Store saved correctly
# get 1 plane
sel = {}
for key, value in  image.im.coords.items():
    sel[key] = value[0]
plane = image.im.sel(sel) 
mean_test = plane.mean().values
logger.debug(f'Plane mean = {mean_test}')
assert mean_test > 0

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


# Process Image    
image.correct_background()
image.focus_projection()
image.register_channels()
if snakemake.params.overlap:
    overlap = int(snakemake.params.overlap)
    direction = snakemake.params.direction
    logger.info(f'Remove {overlap} px {direction} overlap')
    image.remove_overlap(overlap = overlap, direction = direction)
# TODO :: FIX IN pyseq_image 
image.im.name = section_name

# Write Processed Images
save_path = Path(snakemake.output[0]).parents[0]
delayed_store = image.save_zarr(save_path, compute = False)

# Start computation on cluster
logger.info('Processing images')
with performance_report(filename=snakemake.log[1]):
    future_store = client.persist(delayed_store, retries = 10)
    futures = list(future_store.dask.values())
    wait(futures)


    
# Double check no errors
futures_done = [f.done() for f in futures]
if all(futures_done):
    logger.info('Finished processing images')
else:
    logger.info('Error processing images')

    



    ## some dask computation




    

