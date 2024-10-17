from picasso.nn_picasso import PICASSOnn
import xarray as xr
import numpy as np
from pathlib import Path
from dask.distributed import Client, wait, LocalCluster
from math import log
from utils import get_logger
import yaml
from pre import image_analysis as ia
    
section_name = snakemake.params.section

# Start logger
smk_logger = get_logger(section_name, filehandler = snakemake.log[0])
# Open image from zarr store
hs_image = ia.get_HiSeqImages(image_path = snakemake.input[0], logger = smk_logger)
image = hs_image.im
smk_logger.debug(image)

# Make sure only 1 objective step
if 'obj_step' in image.dims:
    if image.obj_step.size > 1:
        mid_step = image.obj_step[image.obj_step.size//2]
        image = image.sel(obj_step = mid_step)

# Get config data
markers_config = snakemake.config.get('markers')
sink_source_ch = snakemake.config.get('unmixing', {610:[558], 740:[687]})
smk_logger.debug(f'sink source channels:: {sink_source_ch}')

# See if background or autofluorescence image exists
AF = ()
for cy, ch_markers in markers_config.items():
    for ch, marker in ch_markers.items():
        if marker in ['background', 'autofluorescence', 'af', 'AF', 'bg', 'BG']:
            AF = (cy,ch)
            smk_logger.info(f'Autofluoresence reference :: cycle {cy} :: channel {ch}')
            



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
        smk_logger.info(f'Unmixing cycle {cy}')
        smk_logger.info(f'Optimizing fluorescence removal from channels {sinks}')
        model = PICASSOnn(cy_mm[cy][0])

        # Fit alpha and background
        for i in model.train_loop(images = image.sel(cycle=cy, channel=all_ch), max_iter=snakemake.params.max_iter):
            pass
        alpha = model.mixing_parameters[0,:,:]
        bg =  model.mixing_parameters[1,:,:]
        picasso_params[cy].update({'alpha':alpha.tolist(), 'background':bg.tolist()})
        smk_logger.info(f'Finished optimizing fluorescence removal from channels {sinks}')
        
# Write picasso parameters to yaml
with open(snakemake.output[0], 'w') as f:
    f.write(yaml.dump(picasso_params))
       
