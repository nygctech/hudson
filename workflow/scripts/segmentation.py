import xarray as xr
import cellpose
from cellpose import core, models, io
from pathlib import Path
import imageio
from utils import get_logger, open_zarr




# Open image from zarr store
image = open_zarr(snakemake.input[0])

# Start logger
smk_logger = get_logger(image.name, filehandler = snakemake.log[0])
smk_logger.debug(image)

# Make sure only 1 objective step
if 'obj_step' in image.dims and 'obj_step' not in snakemake.config.get('segmentation',{}):
    if image.obj_step.size > 1:
        mid_step = image.obj_step[image.obj_step.size//2]
        image = image.sel(obj_step = mid_step)
        smk_logger.debug(image)

# segment
logger = io.logger_setup()
use_GPU = core.use_gpu()
smk_logger.info(f'Using GPU: {use_GPU}')
model = models.CellposeModel(gpu=use_GPU,model_type='TN2')
#model = models.CellposeModel(model_type='TN2')
# Remove once priors steps in pipe
#one_z_plane = image.sel(obj_step = 8498, channel = 558, cycle=1)
arr = image.sel(snakemake.config.get('segmentation'))
channels = [0,0]
masks, flows, styles = model.eval(arr.values, diameter=None, channels=channels, cellprob_threshold= -6, flow_threshold= 1000)
imageio.imwrite(snakemake.output[0],masks)
