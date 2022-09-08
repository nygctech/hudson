from pre import image_analysis as ia
from utils import get_cluster, get_logger
from dask.distributed import Client, wait
from pathlib import Path


# Get section name
section_name = Path(snakemake.input[0]).stem

# Start logger
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])

# Open image
image = ia.get_HiSeqImages(image_path = snakemake.input[0], logname = f'{section_name}.image')

# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
ntiles = int(len(image.im.col)/2048)
min_workers = max(1,2*ntiles)
max_workers = 4*min_workers


# Start computation
with Client(cluster) as client:

    cluster.adapt(minimum = min_workers, maximum=max_workers)
    client.wait_for_workers(min_workers, 60*10)
    
    logger.info(f'Initial name :: {image.im.name}')
    image.correct_background()
    logger.info(f'After correct backround name :: {image.im.name}')
    image.register_channels()
    logger.info(f'After register channels name :: {image.im.name}')
    # TODO :: FIX IN pyseq_image 
    image.im.name = section_name

    # Write Processed Images
    save_path = Path(snakemake.output[0]).parents[0]
    delayed_store = image.save_zarr(save_path, compute = False)
    logger.info('Processing images')
    future_store = client.persist(delayed_store)
    wait(future_store)
   
logger.info('Finished processing images')




    

