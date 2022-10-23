from utils import open_zarr, get_cluster, get_logger
import yaml
import xarray as xr
import dask.array as da
from dask.distributed import Client, wait, performance_report
import dask


# open xarray image from zarr store
image = open_zarr(snakemake.input[0])

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])
logger.info(f'Opened {image.name}')
logger.debug(image.chunks)

# open marker info
markers_config = snakemake.config.get('markers')

# open saved picasso info
with open(snakemake.input[1]) as f:
    picasso_params = yaml.safe_load(f)


# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
nworkers = snakemake.params.tiles
logger.info(f'Scale dask cluster to {nworkers}')
cluster.scale(nworkers)
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

# Write Unmixed Images
unmixed = xr.concat(marker_stack, dim='marker').assign_coords({'marker':marker_name})
#unmixed = unmixed.chunk({'row': len(image.row)})
logger.debug(unmixed)
delayed_store = unmixed.to_dataset().to_zarr(snakemake.output[0], compute = False)
logger.info('Unmixing images')
with performance_report(filename=snakemake.log[1]):
    future_store = client.persist(delayed_store, retries = 10)
    futures = list(future_store.dask.values())
    wait(future_store)
    
futures_done = [f.done() for f in futures]
if all(futures_done):
    logger.info('Finished unmixing images')
else:
    logger.info('Error unmixing images')
