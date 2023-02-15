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
model_type = snakemake.config.get('segmentation',{}).get('model type', 'TN2')
diameter = snakemake.config.get('segmentation',{}).get('diameter', None)
cprob = snakemake.config.get('segmentation',{}).get('cell probability', -6)
fthresh = snakemake.config.get('segmentation',{}).get('flow threshold', 1000)
marker = {'marker':snakemake.config.get('segmentation',{}).get('marker')}
smk_logger.info(f'Segmenting marker')
#TODO Log parameters
smk_logger.info(f'Using model {model_type}')

model = models.CellposeModel(gpu=use_GPU,model_type=model_type)
#model = models.CellposeModel(model_type='TN2')
# Remove once priors steps in pipe
#one_z_plane = image.sel(obj_step = 8498, channel = 558, cycle=1)
sel = snakemake.config.get('segmentation')
for key in sel.copy():
    if key not in image.dims:
        del sel[key]
        
arr = image.sel(sel)
channels = [0,0]
masks, flows, styles = model.eval(arr.values, diameter=diameter, channels=channels, cellprob_threshold= cprob, flow_threshold= fthresh)
imageio.imwrite(snakemake.output[0],masks)
