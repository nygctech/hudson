from pre import image_analysis as ia
from dask import delayed, compute
from dask.distributed import Client, wait, performance_report
import imageio
import numpy as np
from os.path import join
from utils import get_markers

# Start logger
logger = get_logger(image.name, filehandler = snakemake.log[0])

# Open final image
image_path = snakemake.input[0]
im = ia.get_HiSeqImages(image_path)
logger.info(f'Opened {im.im.name}')

col_size = im.im.chunksizes['col'][0]

# Start Cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
logger.debug(cluster.new_worker_spec())
logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
window_size = snakemake.config.get('preprocess', {}).get('focus window', col_size)
ntiles = im.im.col // window_size
logger.info(f'Scale dask cluster to {ntiles}')
cluster.scale(ntiles)
client = Client(cluster)


# Open maske focus map
with np.load(snakemake.input[1]) as npz:
        focus_map = np.ma.MaskedArray(**npz)

# Check markers and make a marker list
nuclei, cells = get_markers()
logger.info(f'Nuclei marker {nuclei}')
logger.info(f'Cell markers {cells}')


window_size = snakemake.config.get('preprocess', {}).get('focus window', col_size)
logger.debug(f'Window size: {window_size}')


# Get threshold parameters for saving coming cell marker images

# Save windows in tile chunk that have tissue
@delayed
def save_windows(image, col_map, c_ind):
    '''Save windows in image coresponding to col_map.

        image = marker and column selected image, ie one chunk
        col_map = masked numpy array, masked values are windows without tissue
        col_ind = column index used for naming images
    '''

    size = image.col.size
    #TODO get marker name
    m = image.marker

    for r, window in enumerate(col_map):
        win_rows = slice(r*size, (r+1)*size)
        if not np.ma.is_masked(window)
            window_im = image.sel(row=win_rows)
            window_name = f'{image.name}_{m}_r{r}c{c_ind}.tiff'
            imageio.imwrite(join(snakemake.output[0],window_name), window_im)

    ## TODO: ADD combined cell marker image
    #if len(cells) > 1:


delayed_windows = []
# loop through markers and then chunks
for m in nuclei + cells:
    for c, col in enumerate(focus_map.T):
        if col.any()
            # focus map is masked where there is no tissue
            tile = im.im.sel(col=slice(c*window_size, (c+1)*window_size), marker = m)
            delayed_windows.append(save_windows(tile, col, c))

logger.info(f'Start saving training images')
with performance_report(filename=snakemake.log[1]):
    client.compute(*delayed_windows)
logger.info(f'Finished saving training images')
