from pre import image_analysis as ia
from dask import delayed, compute
import imageio
import numpy as np
from os.path import join

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])



# Open final image
image_path = snakemake.input[0]
im = ia.get_HiSeqImages(image_path)
logger.info(f'Opened {im.im.name}')

# Start Cluster
# TODO, format cluster and client correctly
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
window_size = snakemake.config.get('preprocess', {}).get('focus window', col_size)
ntiles = im.im.col // window_size
cluster.scale(ntiles)

# Open maske focus map
with np.load(snakemake.input[1]) as npz:
        focus_map = np.ma.MaskedArray(**npz)

nuclei = snakemake.config.get('segmentation', {}).get('nuclei', [])
cells = snakemake.config.get('segmentation', {}).get('cells', [])

# check nuclei is str, if it is convert it to a tolist
# Check markers and make a marker list
if nuclei is a string:
    nuclei = [nuclei]
if cell is a string:
    cell = [cell]
for m in nuclei + cells:
    assert m in im.im.markers.values, f'{m} not listed as a marker'
logger.debug(f'Nuclei marker: {nuclei}')
logger.debug(f'Nuclei marker: {cells}')

col_size = im.im.chunksizes['col'][0]
window_size = snakemake.config.get('preprocess', {}).get('focus window', col_size)
logger.debug(f'Window size: {window_size}')


# Get threshold parameters for saving coming cell marker images

# Save windows in tile chunk that have tissue
@delayed
def save_windows(image, col_mask, col_ind, nuclei=[], cells = []):
    '''Save windows in image coressponding to col_mask of nuclei and cell markers.

        col_mask = masked numpy array, masked values are windows without tissue
        col_ind = column index used for naming images


    '''

    size = image.col.size

    for r, window in enumerate(col_mask):
        if not np.ma.is_masked(window):
            win_rows = slice(r*size, (r+1)*size)
            # loop through nuclei and cell markers
            for m in markers:
                window = image.sel(row=win_rows, marker = m)
                window_name = f'{image.name}_{m}_r{r}c{col_ind}.tiff'
                imageio.imwrite(join(snakemake.output[0],window_name), window)

    ## TODO: ADD combined cell marker image
    #if len(cells) > 1:


delayed_windows = []
# loop through chunks
for c, col in enumerate(focus_map.T):
    if col.any()
        # focus map is masked where there is no tissue
        tile = im.im.sel(col=slice(c*window_size, (c+1)*window_size))
        delayed_windows.append(save_windows(tile, col, c, nuclei = nuclei, cells = cells))

logger.info(f'Start saving training images')
compute(*delayed_windows) # submit to client instead
logger.info(f'Finished saving training images')
