from utils import get_logger, HiSeqImage
import xarray as xr
import cellpose
from cellpose import core, models, io
import imageio
import subprocess
import numpy as np
import yaml
from pathlib import Path
from scipy.signal import find_peaks
from skimage.segmentation import watershed


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
    summary_path = Path(snakemake.output[0]).parents[1] / f"summary_{section_name}.yaml"
    with open(summary_path, 'r') as f:
        summary = yaml.safe_load(f)
        o = summary.get("best_obj_step", None)
        
    if o is not None:
        cytoplasm['obj_step'] = o
    elif image.obj_step.size > 1:
        cytoplasm['obj_step'] = image.obj_step[image.obj_step.size - 1]
    else:
        cytoplasm['obj_step'] = image.obj_step[0]
    
    smk_logger.info(f'Using objective step {cytoplasm["obj_step"]}')


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

# Expand masks for cell+ segments
find_peaks_kws = snakemake.config.get('segmentation').get('expand', {})
if len(find_peaks_kws) > 0: 
    smk_logger.info('Expanding segments')
    smk_logger.info(f"{find_peaks_kws}")
    _im1 = _im1.values
    # Threshold image so we don't label background
    counts, bins = np.histogram(_im1.flatten(), bins=range(0, 4097))
    peaks, props = find_peaks(counts, **find_peaks_kws)
    mask = _im1 > (peaks[0] + 3*props["widths"][0])
    # Expand cell segemnts to cover tissue
    expanded_segments = watershed(-_im1, masks, mask=mask)
    # Save expanded segments
    fname = Path(snakemake.output[0])
    fname = fname.with_name(f"expanded_{section_name}.tiff")
    imageio.imwrite(fname, expanded_segments)

