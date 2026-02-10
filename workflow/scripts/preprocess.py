from pre import image_analysis as ia
from utils import get_cluster, get_logger
from dask.distributed import Client, wait, performance_report
from pathlib import Path


import xarray as xr
import dask.array as da
import numpy as np
from dask_image.ndinterp import affine_transform
from skimage.registration import phase_cross_correlation
import yaml
import imageio
from io import BytesIO

def optimal_obj_step(image):
    im_max = image.max().values
    im_min = image.min().values
    image255 = ((image-im_min)/(im_max - im_min)*255).astype("uint8")
    
    jpeg_size = []
    for o in image255.obj_step.values:
        im = image255.sel(obj_step=o)
        with BytesIO() as f:
            imageio.imwrite(f, im, format='jpeg')
            jpeg_size.append(f.__sizeof__())

    o_ind = np.array(jpeg_size).argmax()
    
    return int(image.obj_step[o_ind].values)


# Get section name
section_name = Path(snakemake.input[0]).stem

# Get Processing Options
focus_projection = snakemake.config.get('preprocess',{}).get('focus projection',False)
registration = snakemake.config.get('preprocess',{}).get('registration',{})

# Start logger
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])
logger.info(f'path:: {snakemake.input[0]}')
logger.info(f'section:: {section_name}')

# Open image
image = ia.get_HiSeqImages(image_path = snakemake.input[0], logname = f'{section_name}.image')
logger.info(f'machine::{image.machine}')
logger.debug(f'{image.im.shape}')

# Check Raw Store saved correctly
# get 1 plane   
sel = {}
for key, value in  image.im.coords.items():
    try:
        sel[key] = value[0]
    except IndexError:
        logger.debug(f'Only 1 {key}')
        pass
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
if focus_projection:
    image.focus_projection()
image.register_channels()

# Register Cycles
if len(registration) > 0:
    cluster.scale(nworkers*2)
    ref_cy = registration.get("reference_cycle", 0)
    ref_ch = registration.get("reference_channel", 610)
    # o = image.im.obj_step[len(image.im.obj_step)//2]
    ref = image.im.sel(cycle=ref_cy, channel=ref_ch)

    # Find optimal obj_step 
    o = optimal_obj_step(ref)
    logger.debug(f"Using obj_step {o} as optimal focus plane")
    ref = ref.sel(obj_step=o)
    summary_path = Path(snakemake.output[0]).parents[1] / f"summary_{section_name}.yaml"
    with open(summary_path, 'r') as f:
        summary = yaml.safe_load(f)
    summary["best_obj_step"] = o
    with open(summary_path, 'w') as f:
        f.write(yaml.dump(summary))
    
    cycles = list(image.im.cycle.values)
    logger.info(f"Registering cycles against cycle {ref_cy} and {ref_ch} nm channel")
    
    def shift_to_affmat(shift):
        aff_matrix = np.identity(3)
        aff_matrix[0,2] = -shift[0]
        aff_matrix[1,2] = -shift[1]
        return aff_matrix
    
    def cycle_affine_transform(im_cy, aff_matrix):
        ch_stack = []
    
        for ch in im_cy.channel:
            o_stack = []
            for o in im_cy.obj_step:
                selection = {"channel":ch, "obj_step":o}
                o_stack.append(affine_transform(im_cy.sel(**selection).data, aff_matrix))
            ch_stack.append(da.stack(o_stack))
                 
        return xr.DataArray(da.stack(ch_stack), dims=im_cy.dims, coords=im_cy.coords)

    cy_stack = []
    for cy in cycles:
        if cy == ref_cy:
            cy_stack.append(image.im.sel(cycle=cy))
        else:
            logger.info(f"Detecting shift in cycle {cy}")
            shift = phase_cross_correlation(ref, image.im.sel(cycle=cy, channel=ref_ch, obj_step=o), upsample_factor=100)
            logger.info(f"Cycle {cy}: shift = {shift[0]}, error = {shift[1]}, phasediff = {shift[2]}")
            aff_matrix = shift_to_affmat(shift[0])
            cy_stack.append(cycle_affine_transform(image.im.sel(cycle=cy), aff_matrix))
    registered = xr.concat(cy_stack, dim="cycle")
    attrs = image.im.attrs
    image.im = registered
    image.im.attrs.update(attrs)
    
# Remove tiling overlap    
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

# save attributes
#with open(save_path.with_suffix('.yaml'), 'w') as f:
#    yaml.dump(image.im.attrs, f)

    

