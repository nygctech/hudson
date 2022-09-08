from utils import open_zarr, get_cluster, get_logger
import yaml
import xarray as xr
import dask.array as da
from dask.distributed import Client, wait


# open xarray image from zarr store
image = open_zarr(snakemake.input[0])

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])
logger.info(f'Opened {image.name}')

# open stain info
stains_config = snakemake.config.get('stains')

# open saved picasso info
with open(snakemake.input[1]) as f:
    picasso_params = yaml.safe_load(f)



# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
ntiles = int(len(image.col)/2048)
min_workers = max(1,2*ntiles)
max_workers = 4*min_workers

dims = ['channel','sink']
with Client(cluster) as client:
    
    cluster.adapt(minimum = min_workers, maximum=max_workers)
    client.wait_for_workers(min_workers, 60*10)
    
    stain_stack = []
    stain_name = []
    for cy, ch_stains in stains_config.items():
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
            ims = image.sel(cycle=cy, channel=all_ch).stack(dims = ['row','col']).T
            ch_stack = []
            for sink_ch in sinks:
                unmixed_im = ((ims - bg.sel(sink=sink_ch)) @ alpha.sel(sink=sink_ch)).unstack().drop_vars('sink')
                ch_stack.append(unmixed_im)
            unmixed_ims = xr.concat(ch_stack, dim='channel').assign_coords({'channel':sinks})

        # Save as markers instead of cycle & channel
        for ch, marker in ch_stains.items():
            if ch in sinks:
                stain_stack.append(unmixed_ims.sel(channel=ch))                     # add unmixed sink
            else:
                stain_stack.append(image.sel(cycle=cy, channel=ch))              # add channel that did not have any spillover
            stain_name.append(marker)

   # Write Unmixed Images
    unmixed = xr.concat(stain_stack, dim='marker').assign_coords({'marker':stain_name})
    delayed_store = unmixed.to_dataset().to_zarr(snakemake.output[0], compute = False)
    logger.info('Unmixing images')
    future_store = client.persist(delayed_store)
    wait(future_store)
    logger.info('Finished unmixing images')
