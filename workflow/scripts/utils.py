def get_cluster(manager='SLURM', **winfo):
    '''Create dask cluster
    
        **Parameters:**
        - manager (str): Name of cluster scheduler
        - winfo: arguments to pass to workers

       **Returns:**
        - dask cluster
    '''

    import re

    assert manager.upper() in ['SLURM']

    if re.search(manager, 'SLURM', re.IGNORECASE):
        from dask_jobqueue import SLURMCluster
        cluster = SLURMCluster(**winfo)

    #Add elif statements to add your specific cluster scheduler

    return cluster


def get_logger(logname = None, filehandler = None):
    '''Creater logger.
    
       **Parameters:**
        - logname (str): Name of logger (optional)
        - filehandler (file object): File object to write log to (optional)

       **Returns:**
        - logging.Logger 
    
    '''
    
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


class HiSeqImage():
    """HiSeqImages

       **Attributes:**
        - im: image as xarray DataArray
        - config: PySeq YAML configuration file from ~/.config/.pyseq2500/machine_settings.yaml
        - machine: Name of machine
        - name: Name of image

    """

    def __init__(self, image_path=None,  im=None, logger=None, **kwargs):
        """The constructor for HiSeq Image Datasets.

           **Parameters:**
            - image_path (path): Path to images stored as zarr, optional
            - im (DataArray): HiSeqImage xarray DataArray, optional
            - logger: Logger for output, optional
            - kwargs: arguments for get_logger

           **Returns:**
            - HiSeqImage object

        """
        
        if logger is None:
            logger = get_logger(**kwargs)
        self.logger = logger

        self.im = None
        self.config  = None
        self.machine = None
        self.name = None
        self.ome_ch_dict = {558: {'excitation_wavelength':532, 'emission_wavelength':558, 'color':'#00429d'},
                            610: {'excitation_wavelength':532, 'emission_wavelength':610, 'color':'#96ffea'},
                            687: {'excitation_wavelength':660, 'emission_wavelength':687, 'color':'#93003a'},
                            740: {'excitation_wavelength':660, 'emission_wavelength':740, 'color':'#ff005e'},
                            }
        self.dim_map = {'marker':'C', 'channel':'C', 'cycle':'T', 'obj_step':'Z', 'row':'Y', 'col':'X'}

        if im is None:
        
            if isinstance(image_path, str):
                image_path = Path(image_path)
            
            #Open zarr image
            im = self.open_zarr(image_path)
                    
            
        if im is not None:
            self.im = im
            self.machine = im.machine
            self.name = im.name
            self.config = get_machine_config(im.machine)
            self.logger.info(f'Opened {im.name}')

    def open_zarr(self, image_path):
        """Create labeled dataset from zarrs.
           
           Xarray wrapped/labeled arrays are stored in self.im list.
    
           **Returns:**
           - list: Names of labeled image datasets 
    
        """
    
    
    
        ome_metadata = zarr.open_group(image_path, mode = 'r').attrs.get('omero',None)
        if ome_metadata is not None:
            # Read zarr written by ome_zarr
            im = read_ome_zarr(image_path, ome_metadata)
        else:
            im = read_xr_zarr(image_path)
    
    
        return im

    def write_ome_zarr(self, dir_path):
        """Write OME zarr from xarray labeled dask data array.

           delayed_write = image.write_ome_zarr('data/')
           dask.compute(*delayed_write)

           **Parameters:**
            - dir_path(str,path): path to directory to store data in

           **Returns:**
           - list of dask delayed objects to write array to disk

        """
        
        if isinstance(dir_path, str):
            dir_path = Path(dir_path)

        # Instrument OME Metadata
        instrument_ = Instrument(microscope = Microscope(**self.config['Microscope']), 
                                 objectives = [Objective(**self.config['Objective'])])
        
        # Channel OME metadata
        channels_ = []
        # print(f'{img.name} dims = {img.im.dims}')
        if 'channel' in self.im.dims:
            for i, ch in enumerate(self.im['channel']):
                # print(f'ID {i}: channel {ch}')
                channels_.append(Channel(name=str(int(ch)), **self.ome_ch_dict[int(ch)]))
        elif 'marker' in self.im.dims:
            for i, ch in enumerate(self.im['marker']):
                # print(f'ID {i}: channel {ch}, cycle {ch.cycle}, em_ch {ch.channel}')
                channels_.append(Channel(name=f'{ch}-Cycle{ch.cycle}', **self.ome_ch_dict[int(ch.channel)]))
        else:
            raise KeyError('Missing Channel Dimension')

        #Dimension Order
        dim_order = ''
        for d in self.im.dims:
            dim_order += self.dim_map.get(d)

        # Pixel OME metadata
        pxs = Pixels(channels = channels_, 
                     dimension_order = 'XYZCT', #Place holder, actually dim order TCZYX
                     size_x = len(self.im.col),
                     size_y = len(self.im.row),
                     size_z = len(self.im.obj_step),
                     size_c = len(self.im.channel),
                     size_t = len(self.im.cycle),
                     tiff_data_blocks = [TiffData(first_z = self.im.obj_step.values[0], first_t = self.im.cycle.values[0])],
                     **self.config['Pixels']
                    )
        
        description =f"""first_cycle = {self.im.cycle[0]},
                        last_cycle = {self.im.cycle[-1]},
                        first_objstep = {self.im.obj_step[0]},
                        last_objstep = {self.im.obj_step[-1]},
                        int_objstep = {self.im.obj_step[1]-self.im.obj_step[0]}
                      """
        ome = OME()
        ome.images = [Image(name = self.name, pixels = pxs, description = description)]
        ome.instruments = [instrument_]
        ome.creator = __name__
        ome_dict = to_dict(ome)

        #Fix UnSerializable Fields
        # Color Object Not JSON Serializable as dictionary
        nch = len(ome_dict['images'][0]['pixels']['channels'])
        for i in range(nch):
            color_val = ome_dict['images'][0]['pixels']['channels'][i]['color'].as_rgb_tuple()
            ome_dict['images'][0]['pixels']['channels'][i]['color'] = color_val
            
        # Pixel Order Object Not JSON Serializable as dictionary
        ome_dict['images'][0]['pixels']['dimension_order'] = dim_order
        ome_dict['images'][0]['pixels']['type'] = ome_dict['images'][0]['pixels']['type'].value
        
        # write the image data
        try:
            zarr_name = dir_path/f'{self.im.name}.zarr' 
            mkdir(zarr_name)
            store = parse_url(zarr_name, mode="w").store
            root = zarr.group(store=store)
            root.attrs["omero"] = ome_dict
            delayed_ = write_image(image=self.im.data, group=root, axes=dim_order.lower(), compute = False,
                                   storage_options = {'chunks':[i[0] for i in self.im.chunks]})
        except Exception as error:
            self.logger.error("Error writing ome_zarr %s", error, exc_info=True)
            # print("Error writing ome_zarr %s", error)
            return False

        return delayed_



def read_xr_zarr(image_path):
    '''Read zarr written by xarray.

       **Parameters**
       image_path (str,path): Path to zarr image data

       **Returns**
       xarray DataArray backed by dask arrays
    
    '''
    
    # Read zarr written by xarray
    im_name = image_path.stem.split('.')[0]
    im = xr.open_zarr(image_path).to_array()
    im = im.squeeze().drop_vars('variable').rename(im_name)

    # read attributes
    attr_path = image_path.with_suffix('.attrs')
    if not attr_path.exists():
        attr_path = attr_path.parents[1] / 'raw_zarr' / attr_path.name

    # add attributes
    try:
        with open(attr_path) as f:
            for line in f:
                items = line.split()
                im.attrs[items[0]] = items[1]
    except:
        print('No attributes read')

    return im
        



def read_ome_zarr(image_path, ome_metadata):
    '''Read zarr written by ome_zarr.

       **Parameters**
       image_path (str,path): Path to zarr image data
       ome_metadata (dict): OME metadata as a dictionary object

       **Returns**
       xarray DataArray backed by dask arrays
    
    '''

    reader = Reader(parse_url(image_path, mode="r"))
    # nodes may include images, labels etc
    nodes = list(reader())
    # first node will be the image pixel data, first item in data list is full res
    im = nodes[0].data[0]
    
    # Read OME Metadata
    im_name = ome_metadata['images'][0]['name']
    meta_dict = {}
    for field in ome_metadata['images'][0]['description'].split(','):
        fn, val = field.split('=')
        meta_dict[fn.strip()] = int(val)
    
    # Map channel, cycle, and obj_step coordinates
    coord_dict = {'C': [c['name'] for c in ome_metadata['images'][0]['pixels']['channels']]}
    if meta_dict.get('first_cycle', False):
        coord_dict['T'] = range(meta_dict['first_cycle'], meta_dict['last_cycle']+1)
    if meta_dict.get('first_objstep', False):
        coord_dict['Z'] = range(meta_dict['first_objstep'], meta_dict['last_objstep']+1, meta_dict['int_objstep'])
    
    dim_map = {'X': 'row', 'Y':'col', 'Z': 'obj_step', 'T': 'cycle', 'C': 'channel'}
    
    dims_ = []; coords_ = {}
    # loop through dimensions ie XYZCT
    for d in ome_metadata['images'][0]['pixels']['dimension_order']: 
        d_name = dim_map[d]
        dims_.append(d_name)
        if d in 'CTZ':
            coords_[d_name] = coord_dict[d]
    
    xr_im = xr.DataArray(data = im, coords = coords_, dims = dims_, name = im_name)
    xr_im.attrs['omero'] = ome_metadata
    xr_im.attrs['machine'] = ome_metadata['instruments'][0]['microscope']['serial_number'].split(',')[0].strip()

    return xr_im


# def open_zarr(image_path):
#     '''Open xarray image from zarr store.'''
    
#     import xarray as xr
#     from pathlib import Path

#     image_path = Path(image_path)
#     im_name = image_path.stem
#     image = xr.open_zarr(image_path).to_array()
#     image = image.squeeze().drop_vars('variable').rename(im_name)
    
#     return image
