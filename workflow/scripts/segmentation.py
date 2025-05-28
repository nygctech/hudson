from utils import get_logger, HiSeqImage
import xarray as xr
import cellpose
from cellpose import core, models, io
import imageio
import subprocess
import numpy as np
import yaml


# Start logger
section_name = snakemake.params.section
smk_logger = get_logger(section_name, filehandler = snakemake.log[0])
smk_logger.debug(section_name)

# Open image from zarr store
hs_image = HiSeqImage(image_path = snakemake.input[0], logger = smk_logger)
image = hs_image.im
smk_logger.debug(image)


# Make sure only 1 objective step
cytoplasm = snakemake.config.get('segmentation').get('cytoplasm')
nuclei = snakemake.config.get('segmentation').get('nuclei', None)
smk_logger.debug(cytoplasm)
smk_logger.debug(nuclei)

for key in cytoplasm.copy():
    if key not in image.dims:
        smk_logger.warning(f'No dimension named {key}')
        del cytoplasm[key]

if 'obj_step' in image.dims and 'obj_step' not in cytoplasm:
    if image.obj_step.size > 1:
        step = image.obj_step[image.obj_step.size - 1]
        cytoplasm['obj_step'] = step
        smk_logger.info(f'Using objective step {step}')       


# segment
logger = io.logger_setup()
use_GPU = core.use_gpu()
res = subprocess.run(['squeue','--format="%.18i %.9P %.30j %.8u %.1T %.10M %.9l %.6D %R"', '--me'], 
                     capture_output=True, text=True)
que = subprocess.run(['grep', f'{section_name}'], input=res.stdout, capture_output=True, text=True)
smk_logger.info(f'{que}')
smk_logger.info(f'Using GPU: {use_GPU}')
seg_args = snakemake.config.get('segmentation')
model_type = seg_args.get('model_type', 'TN2')
diameter = seg_args.get('diameter', 30)
cp_args = {}
cp_args['cellprob_threshold'] = seg_args.get('cell probability', -6)
cp_args['flow_threshold'] = seg_args.get('flow threshold', 1000)


# smk_logger.info(f'Using model {model_type}')
# smk_logger.info(f'diameter = {diameter}')
# smk_logger.info(f'cell probability = {cprob}')
# smk_logger.info(f'flow threshold = {fthresh}')

smk_logger.info(f'Segmenting cytoplasm {cytoplasm}')
_im1 = image.sel(cytoplasm).max('channel')
smk_logger.debug('cytoplasm')
smk_logger.debug(_im1)
if nuclei is None:
    nchan = 1
    im = _im1
    cp_args['channels'] = [0,0]
else:
    nchan = 2
    cp_args['channels'] = [1,2]
    for key in nuclei.copy():
        if key not in image.dims:
            smk_logger.warning(f'No dimension named {key}')
            del nuclei[key]
    if 'obj_step' not in nuclei:
        nuclei['obj_step'] = cytoplasm['obj_step']
    smk_logger.info(f'Segmenting nuclei {nuclei}')
    _im2 = image.sel(nuclei).max('channel')
    smk_logger.debug('nuclei')
    smk_logger.debug(_im2)
    im = xr.concat([_im1, _im2], dim='channel')
    smk_logger.debug('cytoplasm + nuclei')
    smk_logger.debug(im)


model = models.CellposeModel(gpu=use_GPU, model_type=model_type, diam_mean=diameter)
#model = models.CellposeModel(model_type='TN2')
# Remove once priors steps in pipe
#one_z_plane = image.sel(obj_step = 8498, channel = 558, cycle=1)
#sel = snakemake.config.get('segmentation')
# for key in marker.copy():
#     if key not in image.dims:
#         smk_logger.warning(f'No dimension named {key}')
#         del marker[key]
     
# arr = image.sel(marker)
# smk_logger.debug(arr)
# arr = arr.max(dim='channel')
# smk_logger.debug(arr)


cp_args['channel_axis'] = im.dims.index('channel')
smk_logger.info(cp_args)
smk_logger.info('Starting segmentation')
masks, flows, styles = model.eval(im.values, **cp_args)


smk_logger.info('Finished segmentation, calculating metrics')
diam = cellpose.utils.diameters(masks)
std_diam = np.std(diam[1])
mean_diam = np.mean(diam[1])

smk_logger.info('Writing metrics to Summary')
data = {}
data['segmentation'] = {}
data['segmentation']['number_of_cells'] = int(diam[1].shape[0])
data['segmentation']['avergae_cell_size'] = round(float(mean_diam), 4)
data['segmentation']['cell_size_standard_deviation'] = round(float(std_diam), 4)

with open(snakemake.output[1], 'w') as file:
    yaml.dump(data, file)

smk_logger.info('Writing mask')
imageio.imwrite(snakemake.output[0],masks)
