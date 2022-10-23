def get_cluster(manager='SLURM', **winfo):

    import re

    if re.search(manager, 'SLURM', re.IGNORECASE):
        from dask_jobqueue import SLURMCluster as Cluster
    elif re.search(manager, 'HTCondor', re.IGNORECASE):
        from dask_jobqueue import HTCondorCluster as Cluster
    elif re.search(manager, 'LSF', re.IGNORECASE):
        from dask_jobqueue import LSFCluster as Cluster
    elif re.search(manager, 'MOAB', re.IGNORECASE):
        from dask_jobqueue import MOABCluster as Cluster
    elif re.search(manager, 'OAR', re.IGNORECASE):
        from dask_jobqueue import OARCluster as Cluster
    elif re.search(manager, 'PBS', re.IGNORECASE):
        from dask_jobqueue import PBSCluster as Cluster
    elif re.search(manager, 'SGE', re.IGNORECASE):
        from dask_jobqueue import SGECluster as Cluster
    else:
        raise ValueError(f'{manager} not recognized, see https://jobqueue.dask.org/en/latest/api.html')

    cluster = Cluster(**winfo)

    return cluster


def get_logger(logname = None, filehandler = None):
    '''Get logger.'''

    import logging

    if logname is  None:
        logname = __name__

    logger = logging.getLogger(logname)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if filehandler is not None:
        # create file handler which logs even debug messages
        fh = logging.FileHandler(filehandler)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    # . in logname means it's a child logger and don't need to set up console handler
    if '.' not in logname:
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return logger


def open_zarr(image_path):
    '''Open xarray image from zarr store.'''

    import xarray as xr
    from pathlib import Path

    image_path = Path(image_path)
    im_name = image_path.stem
    image = xr.open_zarr(image_path).to_array()
    image = image.squeeze().drop_vars('variable').rename(im_name)

    return image


def get_markers(config):
    nuclei = snakemake.config.get('segmentation', {}).get('nuclei', [])
    cells = snakemake.config.get('segmentation', {}).get('cells', [])
    marker_dict = snakemake.config.get('markers')

    if isinstance(nuclei, str):
        nuclei = list(nuclei)
    if isinstance(cells, str):
        cells = list(cells)

    markers = []
    for cy, chs in marker_dict.items():
        for ch, m in chs.items():
            markers.append(m)

    for m in nuclei:
        assert m in markers, f'nuclei marker, {m}, not listed as a marker'
    for m in cells:
        assert m in markers, f'cell marker, {m}, not listed as a marker'

    return nuclei, cells
