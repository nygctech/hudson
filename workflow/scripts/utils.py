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

def position(pos_):
    """Returns stage position information.

       The center of the image is used to bring the section into focus
       and optimize laser intensities. Image scans of sections start on
       the upper right corner of the section. The section is imaged in
       strips 0.760 mm wide by length of the section long until the entire
       section has been imaged. The box region of interest surrounding the
       section is converted into stage and imaging details to scan the
       entire section.

       =========  ==============================================
         key      description
       =========  ==============================================
       x_center   The xstage center position of the section.
       y_center   The ystage center position of the section.
       x_initial  Initial xstage position to scan the section.
       y_initial  Initial ystage position to scan the section.
       x_final    Last xstage position of the section scan
       y_final    Last ystage position of the section scan
       n_tiles    Number of tiles to scan the entire section.
       n_frames   Number of frames to scan the entire section.
       =========  ==============================================

       **Parameters:**
        - pos(str): AorB: LLx, LLy, URx, URy 
                    where AorB = A or B, LL=Lower Left, UR=Upper Right,  corner using the slide ruler


        For Example: 'A: 14, 20.5, 11, 18'
 
          Left and UR=Upper Right corner using the slide ruler.

       **Returns:**
        - dict: Dictionary of stage positioning and imaging details to scan
          the entire section. See table above for details.
    """

    AorB, box = pos_.split(':')
    AorB = AorB.strip()
    box = [float(_) for _ in box.split(',')]
    
    pos = {}

    LLx = box[0]
    LLy = box[1]
    URx = box[2]
    URy = box[3]

    ####TODO read in from machine_settings.yaml
    tile_width = 0.769                                                 #mm
    resolution = 0.375    # um
    x_spum = 0.4096     #steps per um
    y_spum = 100     # steps per um
    fc_origin = {'A':[17571,-180000], 'B':[43310,-180000]}
    bundle_height = 128
    overlap = 0; overlap_dir = 'left'
    #####
    
    # Number of scans
    dx = tile_width-resolution*overlap/1000                  # x stage delta in in mm
    n_tiles = ceil((LLx - URx)/dx)
    pos['n_tiles'] = n_tiles

    # X center of scan
    x_center = fc_origin[AorB][0]
    x_center -= LLx*1000*x_spum
    x_center += (LLx-URx)*1000/2*x_spum
    x_center = int(x_center)

    # initial X of scan
    x_initial = n_tiles*dx*1000/2                                           #1/2 fov width in microns
    if overlap_dir == 'left':
        x_initial -= resolution*overlap                           #Move stage to compensate for discarded initial px
    x_initial = int(x_center - x_initial*x_spum)
    pos['x_initial'] = x_initial

    # initial Y of scan
    y_initial = int(fc_origin[AorB][1] + LLy*1000*y_spum)
    pos['y_initial'] = y_initial

    # Y center of scan
    y_length = (LLy - URy)*1000
    y_center = y_initial - y_length/2*y_spum
    y_center = int(y_center)

    # Number of frames
    n_frames = y_length/bundle_height/resolution
    pos['n_frames'] = ceil(n_frames + 10)

    # Adjust x and y center so focus will image (32 frames, 128 bundle) in center of section
    x_center -= int(tile_width*1000*x_spum/2)
    pos['x_center'] = x_center
    y_center += int(32*bundle_height/2*resolution*y_spum)
    pos['y_center'] = y_center

    # Calculate final x & y stage positions of scan
    pos['y_final'] = int(y_initial - y_length*y_spum)
    pos['x_final'] = int(x_initial +(LLx - URx)*1000*x_spum)
    pos['obj_pos'] = None

    return pos
