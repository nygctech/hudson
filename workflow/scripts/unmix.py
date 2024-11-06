
from utils import get_cluster, get_logger
import yaml
import xarray as xr
import dask.array as da
from dask.distributed import Client, wait, performance_report
import dask
from pathlib import Path
from math import floor
from pre import image_analysis as ia
from os import makedirs 


section_name = snakemake.params.section
logger = get_logger(section_name, filehandler = snakemake.log[0])

# Open image from zarr store
hs_image = ia.get_HiSeqImages(image_path = snakemake.input[0], logger = logger)
image = hs_image.im

# Start logger
# logger.info(f'Opened {image.name}')
logger.debug(image.chunks)

# open marker info
markers_config = snakemake.config.get('markers')

# open saved picasso info
with open(snakemake.input[1]) as f:
    picasso_params = yaml.safe_load(f)
    
# open section summary
with open(snakemake.input[2]) as f:
    section_summary = yaml.safe_load(f)


# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
logger.info(f'Dask worker settings: {winfo}')
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
nworkers = section_summary['section_information'].get('tiles', 2)
logger.info(f'Scale dask cluster to {nworkers}')
cluster.scale(nworkers)
logger.info(f'Dask cluster info {cluster}')
client = Client(cluster)

# Save channel/cycle as markers and unmix
dims = ['channel','sink']
marker_stack = []
marker_name = []
for cy, ch_markers in markers_config.items():
    params = picasso_params.get(cy)
    sinks = params.get('sinks')

    if len(sinks) > 0:
        logger.info(f'Cycle {cy} :: sinks {sinks}') 
        all_ch = params.get('channels')
        logger.info(f'Cycle {cy} :: channels {all_ch}') 
        coords = {'channel': all_ch, 'sink':sinks}
        alpha = xr.DataArray(da.from_array(params.get('alpha')), name = 'alpha', dims = dims, coords = coords)
        bg = xr.DataArray(da.from_array(params.get('background')), name = 'background', dims = dims, coords = coords)

        # Unmix Images
#         with dask.config.set(**{'array.slicing.split_large_chunks': False}):
#             ims = image.sel(cycle=cy, channel=all_ch).stack(px = ('row','col'))
        ims = image.sel(cycle=cy, channel=all_ch)
        ch_stack = []
        for sink_ch in sinks:
            unmixed_im = (ims - bg.sel(sink=sink_ch)).clip(min=0).dot(alpha.sel(sink=sink_ch))
            unmixed_im = unmixed_im.clip(min=0).astype('uint16')
#             unmixed_im = unmixed_im.unstack().drop_vars('sink')
            unmixed_im = unmixed_im.drop_vars('sink')
            ch_stack.append(unmixed_im)
        unmixed_ims = xr.concat(ch_stack, dim='channel').assign_coords({'channel':sinks})

    # Save as markers instead of cycle & channel
    for ch, marker in ch_markers.items():
        if ch in sinks:
            marker_stack.append(unmixed_ims.sel(channel=ch))                     # add unmixed sink
        else:
            marker_stack.append(image.sel(cycle=cy, channel=ch))              # add channel that did not have any spillover
        marker_name.append(marker)


# Make xarray unmixed image
unmixed = xr.concat(marker_stack, dim='marker').assign_coords({'marker':marker_name})
unmixed = unmixed.sel(row=range(64,len(unmixed.row)))
logger.debug(unmixed)
# Make HiSeqImage
unmixed = ia.HiSeqImages(im=unmixed, machine=hs_image.machine, files=[snakemake.input[0]], logger=logger)
save_path = Path(snakemake.output[0]).parent
# Write Unmixed Images
delayed_store = unmixed.save_ome_zarr(save_path)
logger.info('Unmixing images')
with performance_report(filename=snakemake.log[1]):
    logger.debug(delayed_store)
    future_store = client.persist(delayed_store, retries = 10)
    logger.debug(future_store)
    futures = [list(f.dask.values())[0] for f in future_store]
    logger.debug(futures)
    wait(futures)
    logger.debug(futures)
    
# Double check no errors
futures_done = [f.done() for f in futures]
logger.debug(futures_done)
if all(futures_done):
    logger.info('Finished unmixing images')
else:
    logger.info('Error unmixing images')

# Write preview images
downscale = floor(section_summary['section_information']['planesizeMB']/25)
makedirs(snakemake.output[1], exist_ok=True)
unmixed.preview_jpeg(image_path=snakemake.output[1], downscale=downscale)

