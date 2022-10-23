from picasso.nn_picasso import PICASSOnn
import xarray as xr
import numpy as np
from pathlib import Path
from dask.distributed import Client, wait, LocalCluster
from math import log
from utils import get_logger, open_zarr
import yaml
    
# Open image from zarr store
image = open_zarr(snakemake.input[0])

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])
logger.info(f'Opened {image.name}')

# Make sure only 1 objective step
if 'obj_step' in image.dims:
    if image.obj_step.size > 1:
        mid_step = image.obj_step[image.obj_step.size//2]
        image = image.sel(obj_step = mid_step)

# Get config data
markers_config = snakemake.config.get('markers')
sink_source_ch = snakemake.config.get('unmixing', {610:[558], 740:[687]})
logger.debug(f'sink source channels:: {sink_source_ch}')

# See if background or autofluorescence image exists
AF = ()
for cy, ch_markers in markers_config.items():
    for ch, marker in ch_markers.items():
        if marker in ['background', 'autofluorescence', 'af', 'AF', 'bg', 'BG']:
            AF = (cy,ch)
            logger.info(f'Autofluoresence reference :: cycle {cy} :: channel {ch}')
            



# Find which channels to unmix
cy_mm = {}
picasso_params = {}
for cy, ch_markers in markers_config.items():
    si_so = {}                                                                  # sink channel : [source channels]
    all_ch = []                                                                 # all channels that need to be used for unmixing

    # Get dictionary of sink channels that need source channels removed
    for ch in ch_markers.keys():                                                 # loop over channels used in cycle
        for source_ch in sink_source_ch.get(ch,[]):                             # loop over source channels if a channel has spillover fluorescence
            if source_ch in ch_markers.keys():                                   # add source channel to sink_source dict (si_so)
                si_so.setdefault(ch,[]).append(source_ch)
                all_ch.append(source_ch)
    all_ch += si_so.keys()
    all_ch = np.unique(all_ch)

    # Making mixing matrix for cycle
    mm = np.zeros((len(all_ch), len(si_so)))
    for c, (si, so_) in enumerate(si_so.items()):
        si_ind = np.where(all_ch == si)
        mm[si_ind, c] = 1
        for so in so_:
            so_ind = np.where(all_ch == so)
            mm[so_ind, c] = -1

    sinks = list(si_so.keys())
    cy_mm[cy] = [mm, all_ch, sinks]                          # {cycle: [mixing matrix, all images, sinks]}
    picasso_params[cy] = {'sinks': sinks, 'channels':all_ch.tolist()}
    
# Unmix images and stack by marker instead of cycle > channel
marker_stack = []
marker_name = []
dims = ['channel','sink']
for cy , ch_markers in markers_config.items():
    mm, all_ch, sinks = cy_mm.get(cy)

    if len(sinks) > 0:
        logger.info(f'Unmixing cycle {cy}')
        logger.info(f'Optimizing fluorescence removal from channels {sinks}')
        model = PICASSOnn(cy_mm[cy][0])

#         # DELETEME
#         test = image.sel(cycle=cy, channel=all_ch)
#         max_px = test.max().compute()
#         min_px = test.min().compute()
#         M = 1                                                                   # Max MINE model value
#         d = max_px - min_px                                                     # dimension of parameter space
#         K = 1                                                                   # Max MINE parameter value
#         L = 10                                                                  # Lipshitz constant, no idea just overestimate
#         e = 0.1                                                            # accuracy
#         c = 0.9                                                          # confidence

#         print(max_px, min_px)
#         print((2*M**2*(d*log(16*K*L*d**(0.5)/e) + 2*d*M + log(2/c)))/(e**2))

        # Fit alpha and background
        for i in model.train_loop(images = image.sel(cycle=cy, channel=all_ch), max_iter=snakemake.params.max_iter):
            pass
#         logger.info('Mixing params:')
#         coords = {'channel':all_ch, 'sink':sinks}
#         alpha = xr.DataArray(model.mixing_parameters[0,:,:], name = 'alpha', dims = dims, coords = coords)
#         bg = xr.DataArray(model.mixing_parameters[1,:,:], name = 'background', dims = dims, coords = coords)
        alpha = model.mixing_parameters[0,:,:]
        bg =  model.mixing_parameters[1,:,:]
        picasso_params[cy].update({'alpha':alpha.tolist(), 'background':bg.tolist()})
        logger.info(f'Finished optimizing fluorescence removal from channels {sinks}')
        
# Write picasso parameters to yaml
with open(snakemake.output[0], 'w') as f:
    f.write(yaml.dump(picasso_params))

#         # Unmix Images
#         ims = image.sel(cycle=cy, channel=all_ch).stack(dims = ['row','col']).T
#         ch_stack = []
#         for sink_ch in sinks:
#             unmixed_im = ((ims - bg.sel(sink=sink_ch)) @ alpha.sel(sink=sink_ch)).unstack().drop_vars('sink')
#             ch_stack.append(unmixed_im)
#         unmixed_ims = xr.concat(ch_stack, dim='channel').assign_coords({'channel':sinks})

#     # Save as markers instead of cycle & channel
#     for ch, marker in ch_stains.items():
#         if ch in sinks:
#             stain_stack.append(unmixed_ims.sel(channel=ch))                     # add unmixed sink
#         else:
#             stain_stack.append(image.sel(cycle=cy, channel=ch))              # add channel that did not have any spillover
#         stain_name.append(marker)
        
# unmixed = xr.concat(stain_stack, dim='marker').assign_coords({'marker':stain_name})  
# unmixed.to_dataset().to_zarr(snakemake.output[0], compute = True)
   
# # Save images  

# cluster = LocalCluster()
# logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
# cluster.scale(snakemake.resources.cores)

            
# zarr_path = snakemake.output[0]  
# delayed_store = unmixed.to_dataset().to_zarr(zarr_path, compute = False)
# logger.info('Writing Unmixed Images')
# future_store = client.persist(delayed_store)
# wait(future_store)
# logger.info('Finished Writing Unmixed Images')

# client.close()
       
