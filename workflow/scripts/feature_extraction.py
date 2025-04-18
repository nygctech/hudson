import zarr
import xarray as xr
import numpy as np
import torch
import joblib
from skimage.io import imread
import time
import pickle
from os.path import exists, join
from joblib import Parallel, delayed
from joblib import parallel_backend
from utils import get_cluster, get_logger, HiSeqImage
from skimage.measure import regionprops_table
import anndata as ad
import mudata as md
import pandas as pd
import squidpy as sq




### Morphological features computed when making muon object !!!

section_name = snakemake.params.section

# Start logger
smk_logger = get_logger(section_name, filehandler = snakemake.log[0])

# Open image from zarr store
hs_image = HiSeqImage(image_path = snakemake.input[0], logger = smk_logger)
image = hs_image.im
smk_logger.debug(image)



# Open instance labels
labels = imread(snakemake.input[1])
smk_logger.info(f'Opened {image.name} labels')


##########################################################################################
## Morphological Features 
##########################################################################################



# Get Morphological Features
smk_logger.info(f'Measuring morphological features')
features = ('area','area_bbox','area_convex','area_filled','axis_major_length','axis_minor_length',
            'eccentricity', 'equivalent_diameter_area','euler_number','extent','feret_diameter_max',
            'orientation','perimeter','perimeter_crofton','solidity', 'label', 'centroid')
morph_table = regionprops_table(labels, properties=features)

# Organize Morphological Features
coords = np.vstack([morph_table['centroid-0'], morph_table['centroid-1']]).T
ind = morph_table['label']

del morph_table['centroid-0']
del morph_table['centroid-1']
del morph_table['label']


# Make Morphological AnnData object
morph_df = pd.DataFrame(data = morph_table, index = ind, dtype = np.single)
morph_ad = ad.AnnData(X = morph_df, obsm = {'spatial':coords})

# Add graph
radius = None # get from config file
if radius is None:
    radius = (morph_table['feret_diameter_max'].mean()+morph_table['feret_diameter_max'].std())*3/2
smk_logger.info(f'Used {radius} px radius for delaunay graph')
sq.gr.spatial_neighbors(morph_ad, n_neighs = True, delaunay=True, radius=(0,radius), coord_type = 'generic')
# h5 can't save radius parameter as tuple, so save as float
morph_ad.uns['spatial_neighbors']['params']['radius'] = radius

feat_dict = {'morphological': morph_ad}

##########################################################################################
## Protein Features 
##########################################################################################

# Max project if Z dimension still present
if 'obj_step' in image.dims:
    smk_logger.info(f'Projecting Z max')
    image = image.max(dim = 'obj_step')
    
# Get markers
marker_list = list(image.channel.values)
for m in snakemake.config.get('feature_extraction',{}).get('exclude', []):
    marker_list.remove(m)
msg = 'Intensity features: '
for m in marker_list:
    msg += f' {m}'
smk_logger.info(msg)


##TODO: measure texture / intensity quartiles
mean_intensity_per_marker = {}
for m in image.channel.values:
    smk_logger.info(f'Measuring {m}')
    props = regionprops_table(labels, intensity_image = image.sel(channel = m).values, 
                                  properties = ('intensity_mean',))
    mean_intensity_per_marker.update({m:props['intensity_mean']})

# Make Protein object
prot_df = pd.DataFrame(data = mean_intensity_per_marker, index = ind, dtype = np.single)
prot_ad = ad.AnnData(X = prot_df)

# Save intensity histograms 
protein_names = np.array(prot_ad.var.index)
protein = prot_ad.X

df = pd.DataFrame(columns=['Protein Name', 'Counts', 'Bins'])
for i, pro in enumerate(protein_names):
    counts, bins= np.histogram(protein[:,i], bins=20)
    df.loc[i] = [pro, counts, bins]

df.to_csv(snakemake.output[1], index=False)

feat_dict['protein'] = prot_ad        
##########################################################################################
## Imagenet Features 
##########################################################################################

# Check config for imagenet options
extract_imagenet = snakemake.config.get('feature_extraction',{}).get('imagenet', {})
color_dict = {}
try:
    for k in extract_imagenet.keys():
        if k.upper() in 'RGB':
            color_dict[k.upper()] = extract_imagenet[k]
        else:
            smk_logger.info(f'Could not map {extract_imagenet[k]} to RGB')
            smk_logger.info(f'Assign marker to color, R: marker')
except Exception as e:
    smk_logger.info(f'imagenet error: {e}')
    

# Run imagenet
# TODO: run imagenet with 1 or 2 features, or > 3 (but probably have to optimize more or extend resources)
if len(color_dict.keys()) == 3:
    

    import torch
    from timm.data import resolve_data_config
    from timm.data.transforms_factory import create_transform
    import timm
    from skimage.measure import regionprops
    from dask import delayed
    import dask.array as da
    import dask.dataframe as df
    from dask.distributed import Client, wait, performance_report
    from torchvision.transforms.functional import to_pil_image
    from math import ceil
    from pathlib import Path
    import pandas as pd
    
    
    # TODO: Create options in snakmake config, or pull from machine config, or embed zarr attrs
    #dim = 'marker'
    #px_min = 0; px_max= 4095;
    
    # Confirm Markers
    markers_ = []
    for k in 'RGB':
        if k in color_dict.keys():
            smk_logger.info(f'{color_dict[k]} assigned to {k} channel')
            markers_.append(color_dict[k])
    
    ############### Loading (ViT) model from timm package ##############
    smk_logger.info("initializing imagenet model...")
    model = timm.create_model('vit_base_patch16_224_miil.in21k', pretrained=True)
    model.eval()
    config = resolve_data_config({}, model=model)
    transform = create_transform(**config)
    smk_logger.info("imagenet model initialized")
    
    # Start dask cluster
    # specify default worker options in ~/.config/dask/jobqueue.yaml
    winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
    cluster = get_cluster(**winfo)
    smk_logger.debug(cluster.new_worker_spec())
    smk_logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
    ntiles = image.col.size//2048
    nworkers = max(2,ntiles*2*2)
    smk_logger.info(f'Scale dask cluster to {nworkers}')
    cluster.scale(nworkers)
    client = Client(cluster)
    client.wait_for_workers(ceil(nworkers/4))
    
    # Get label bbox and image_filled
    # Might want to do this with regionprops_table
    props = regionprops(labels)
    
    #Run cell through imagenet 
    @delayed
    def get_logits(im, mask, transform, model, px_min = 0, px_max = 4095):

        mask = np.array([mask]*3)

        im = ((im - px_min)/(px_max-px_min)*255).astype('uint8').values
        im[~mask] = 0
        pil = to_pil_image(torch.tensor(im))
        tensor = (transform(pil)).unsqueeze(0)


        with torch.no_grad():
            logits = model(tensor).to('cpu').detach().numpy()

        return logits


    # Crop out image of cell
    def gen_cells(props, im):
        n_cells = len(props)

        n = 0

        while n < n_cells:
            # get cell mask / bounding box
            rmin, cmin, rmax, cmax = props[n].bbox
            cell = im.sel({'row':slice(rmin, rmax), 'col':slice(cmin, cmax)})


            yield props[n], cell

            n += 1

            
    # send model/transform to dask workers
    dask_model = client.scatter(model, broadcast=True)
    dask_transform = client.scatter(transform, broadcast=True)

    # Loop through cells and compute imagenet features
    logit_stack = []
    for p, c in gen_cells(props, image.sel(channel = markers_)):
        logits = get_logits(c, p.image_filled, dask_transform, dask_model)
        logit_stack.append(da.from_delayed(logits, shape = (1,11221), dtype = np.single))    #Shape is only for 1 set of markers
    imagenet_ = da.concatenate(logit_stack).rechunk() 


    # delayed_store = features.to_zarr(Path(snakemake.output[1]), compute = False)
                                     
    # Write imagenet features to file and log cluster performance
    smk_logger.info(f'Computing imagenet features')
    cluster_report = Path(snakemake.log[0]).with_name(f'features_{image.name}.html')
    with performance_report(filename=cluster_report):
        imagenet = imagenet_.compute()
        
    imagenet_df = pd.DataFrame(data = imagenet, index = ind, columns = [f'{i:05d}' for i in range(11221)], dtype = np.single)
    feat_dict['imagenet'] = ad.AnnData(X = imagenet_df)


# Write features to file
smk_logger.info(f'Writing features')
mdata = md.MuData(feat_dict)
mdata.write(snakemake.output[0])
smk_logger.info('Completed writing features')
