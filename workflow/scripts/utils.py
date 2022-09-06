def get_cluster(manager='SLURM', **winfo):

    import re

    assert manager.upper() in ['SLURM']

    if re.search(manager, 'SLURM', re.IGNORECASE):
        from dask_jobqueue import SLURMCluster
        cluster = SLURMCluster(**winfo)

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
