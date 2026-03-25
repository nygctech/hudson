#from snakemake.script import snakemake
import numpy as np
import torch
from skimage.io import imread
from joblib import delayed
from utils import get_cluster, get_logger, HiSeqImage
from skimage.measure import regionprops_table, regionprops
from skimage.feature import local_binary_pattern
import anndata as ad
import mudata as md
import pandas as pd
import squidpy as sq
from pathlib import Path




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

# Get config options
config = snakemake.config.get('feature_extraction',{})


##########################################################################################
## Morphological Features 
##########################################################################################



# Get Morphological Features
smk_logger.info('Measuring morphological features')
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
morph_var = pd.DataFrame({"feature": ["morphology"]*len(morph_table)}, index=morph_table.keys())
morph_ad = ad.AnnData(X=morph_df, var=morph_var, obsm={'spatial':coords})

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
    smk_logger.info('Projecting Z max')
    image = image.max(dim = 'obj_step')
    
# Get markers
marker_list = list(image.channel.values)
for m in config.get('exclude', []):
    marker_list.remove(m)
msg = 'Intensity features: '
for m in marker_list:
    msg += f' {m}'
smk_logger.info(msg)

# Get Texture markers and local binary pattern (lbp) options
texture_markers = config.get("texture",{}).get("markers", [])
if len(texture_markers) == 0:
    texture_markers = marker_list
lbp_radius = config.get("texture",{}).get("radius", 1)
lbp_points = 8*lbp_radius
lbp_bins = np.arange(0, lbp_points + 3)
def lbp_hist(mask_im, intensity_im):
    cell = intensity_im[mask_im]
    (hist, _) = np.histogram(cell, bins=lbp_bins, density=True)
    return hist
def lbp_hist2(props):
    hists = []
    for prop in props:
        cell = prop.image_intensity[prop.image]
        hists.append(np.histogram(cell, lbp_bins, density=True)[0])
    return hists


mean_intensity_per_marker = {}
texture_feats = []
protein_var = []
for m in marker_list:
    smk_logger.info(f'Measuring {m}')
    marker_im = image.sel(channel = m).values
    props = regionprops_table(labels, intensity_image = marker_im, 
                                  properties = ('intensity_mean',))
    mean_intensity_per_marker.update({m:props['intensity_mean']})
    protein_var.append("protein")

prot_df = pd.DataFrame(data = mean_intensity_per_marker, index = ind, dtype = np.single)

# Measure features inside expanded cell segments
if len(snakemake.config.get('segmentation').get('expand', {})) > 0: 
    fname = Path(snakemake.input[1])
    fname = fname.with_name(f"expanded_{section_name}.tiff")
    labels = imread(fname)
    expanded_mean_intensity_per_marker = {}
    for m in marker_list:
        smk_logger.info(f'Measuring {m} in expanded segment')
        marker_im = image.sel(channel=m).values
        props = regionprops_table(labels, intensity_image = marker_im, properties = ('intensity_mean',"label"))
        expanded_mean_intensity_per_marker.update({f"expanded_{m}":props['intensity_mean']})
        protein_var.append("expanded")

        # Measure texture in expanded segment
        if m in texture_markers:
            smk_logger.info(f'Computing local binary pattern of {m}')
            lbp_im = local_binary_pattern(marker_im, lbp_points, lbp_radius, method='uniform')
            smk_logger.info(f'Calculating {m} texture in expanded segments')
            #props = regionprops_table(labels, intensity_image = lbp_im, extra_properties = (lbp_hist,))
            lbp_props = regionprops(labels, intensity_image = lbp_im)
            hists = lbp_hist2(lbp_props)
            col = [f'{m}_t{i+1:02d}' for i in range(lbp_points+2)]
            texture_feats.append(pd.DataFrame(hists, index=props["label"], columns=col))
            protein_var += ["texture"]*(lbp_points+2)
             
    ind = props["label"]
    expanded_prot_df = pd.DataFrame(data = expanded_mean_intensity_per_marker, index = ind, dtype = np.single)
    texture_prot_df = pd.concat(texture_feats, axis=1)
    prot_df = prot_df.join([expanded_prot_df, texture_prot_df]).fillna(0)

protein_var = pd.DataFrame({"feature": protein_var}, index=prot_df.columns)
prot_ad = ad.AnnData(X = prot_df, var=protein_var)
feat_dict['protein'] = prot_ad  

# Save intensity histograms 
protein_names = np.array(prot_ad.var.index)
protein = prot_ad[:, prot_ad.var["feature"] == "protein"].X
df = pd.DataFrame(index = range(0, 4096-16, 16))
for i, m in enumerate(marker_list):
    counts, bins= np.histogram(protein[:,i], bins=range(0, 4096, 16))
    df[m] = counts
df.to_csv(snakemake.output[1])
      
##########################################################################################
## Imagenet Features 
##########################################################################################

# Check config for imagenet options
imagenet_cfg = config.get('imagenet', {})
color_dict = {}
try:
    marker_colors = imagenet_cfg.get("color", {})
    for k in marker_colors.keys():
        if k.upper() in 'RGB':
            color_dict[k.upper()] = marker_colors[k]
        else:
            smk_logger.info(f'Could not map {marker_colors[k]} to RGB')
            smk_logger.info('Assign marker to color, R: marker')
except Exception as e:
    smk_logger.info(f'imagenet error: {e}')
    

# Run imagenet
# TODO: run imagenet with 1 or 2 features, or > 3 (but probably have to optimize more or extend resources)
if len(color_dict.keys()) == 3:
    

    import torch
    from torchvision.transforms import v2
    import timm
    from skimage.measure import regionprops
    from dask import delayed
    import dask.array as da
    import dask.dataframe as df
    from dask.distributed import Client, performance_report
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
    model_name = imagenet_cfg.get("model", "vit_base_patch16_224_miil.in21k")
    smk_logger.info(f"initializing {model_name} model...")
    model = timm.create_model(model_name, 
                              pretrained=True,
                              num_classes=0) #Removes classifier layer
    model.eval()

    # Transform to resize images and rescale pixel values
    def get_transform(config):
        interp_mode = v2.InterpolationMode(config["interpolation"])
        mean = config["mean"]
        std = config["std"]
        input_size = config["input_size"]
        
        transform = v2.Compose([
            v2.ToImage(),
            v2.Resize(size=input_size[-1], interpolation=interp_mode, antialias=True),
            v2.CenterCrop(size=input_size[1:]),
            v2.ToDtype(torch.float32, scale=True),
            v2.Normalize(mean=mean, std=std)
        ])
        
        return transform
    transform = get_transform(model.pretrained_cfg)
    
    smk_logger.info(f"{model_name} initialized")
    
    # Start dask cluster
    # specify default worker options in ~/.config/dask/jobqueue.yaml
    winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
    # For inference, better to use 1 core with more memory
    rescale_workers = winfo.get("cores", 1)
    memory = winfo.get("memory", "16G")
    winfo["cores"] = 1
    winfo["memory"] = f"{int(memory.split("G")[0])*rescale_workers}G"
    cluster = get_cluster(**winfo)
    smk_logger.debug(cluster.new_worker_spec())
    smk_logger.info(f'cluster dashboard link:: {cluster.dashboard_link}')
    ntiles = image.col.size//2048
    nworkers = max(2*rescale_workers,ntiles*2*2)
    smk_logger.info(f'Scale dask cluster to {nworkers}')
    cluster.scale(nworkers)
    client = Client(cluster)
    client.wait_for_workers(ceil(nworkers/4))

    # send model/transform to dask workers
    smk_logger.info('Scatter model')
    dask_model = client.scatter(model, broadcast=True)
    dask_transform = client.scatter(transform, broadcast=True)

    
    #Run cells through ViT 
    @delayed
    def embed_batch(cell_ims, transform, model):
    
        batch_results = []
        for im in cell_ims:
            tensor = transform(im.values)
            batch_results.append(tensor)
        batch_tensor = torch.stack(batch_results)

        with torch.no_grad():
            embeddings = model(batch_tensor).to('cpu').detach().numpy()
            
        return embeddings

    # Crop out image of cell
    def gen_cells(labels, image):
        props = regionprops(labels, image)
        for prop in props:
            yield prop.image_intensity

        #n_cells = len(props)
        #n = 0
        #while n < n_cells:
        #    # get cell mask / bounding box
        #    rmin, cmin, rmax, cmax = props[n].bbox
        #    cell = im.sel({'row':slice(rmin, rmax), 'col':slice(cmin, cmax)})
        #    
        #    yield cell
        #
        #    n += 1

    # Group cells into batches
    def batch_generator(iterable, size):
        batch = []
        for item in iterable:
            batch.append(item)
            if len(batch) == size:
                yield batch
                batch = []
        if batch:
            yield batch

    # Normalize images against expression mean and rescale values between 0 and 1
    def normalize_image(image, markers, protein_ad):

        smk_logger.info('Normalizing RGB image')
    
        image = image.sel(channel=markers).transpose("row", "col", "channel")
        
        if f"expanded_{markers[0]}" in protein_ad.var_names:
            _markers = [f"expanded_{m}" for m in markers_] 
        else:
            _markers = markers
    
        # Mean expression of markers in segments
        mean_exp = prot_ad[:, _markers].X.mean(axis=0)
        smk_logger.info(f'Mean expression in segments {mean_exp}')
        # Max pixel value of markers
        max_px = image.max(dim=["row", "col"]).compute()
        smk_logger.info(f'Max px in image {max_px}')
        max_px = max_px/mean_exp
        # Min expression of markers in segments (~tissue background)
        min_exp = prot_ad[:, _markers].X.min(axis=0)
        smk_logger.info(f'Min expression in segments {min_exp}')
        min_exp = min_exp/mean_exp
    
        # Normalize images against mean expression
        norm_image = image / mean_exp
    
        # Rescale image values betweeen 0 & 1
        # Don't use max expression of markers in segments,
        # it significantly clips the range of pixel values
        scale_image = (norm_image-min_exp)/(max_px - min_exp)
        scale_image = scale_image.clip(min=0.0, max=1.0)
    
        return scale_image

    # Loop through batches of cells and compute imagenet features
    batch_size = imagenet_cfg.get('batch_size', 32)
    smk_logger.info(f'Batching cells: batch size = {batch_size}')
    embedding_stack = []
    norm_image = normalize_image(image, markers_, prot_ad)
    cell_gen = gen_cells(labels, norm_image)
    for batch in batch_generator(cell_gen, batch_size):
        embeddings = embed_batch(batch, dask_transform, dask_model)
        embedding_stack.append(da.from_delayed(embeddings, shape = (len(batch),768), dtype = np.single)) 
    _imagenet = da.concatenate(embedding_stack, axis=0)

    # delayed_store = features.to_zarr(Path(snakemake.output[1]), compute = False)
                                     
    # Write imagenet features to file and log cluster performance
    smk_logger.info('Computing ViT embeddings')
    cluster_report = Path(snakemake.log[0]).with_name(f'features_{image.name}.html')
    with performance_report(filename=cluster_report):
        imagenet = _imagenet.compute()

    imagenet_col = [f'{i:03d}' for i in range(imagenet.shape[1])]
    imagenet_df = pd.DataFrame(data = imagenet, index = ind, columns = imagenet_col, dtype = np.single)
    imagenet_var = pd.DataFrame({"feature":["imagenet"]*len(imagenet_col)}, index=imagenet_col) 
    feat_dict['imagenet'] = ad.AnnData(X = imagenet_df, var=imagenet_var)

# Write features to file
smk_logger.info('Writing features')
mdata = md.MuData(feat_dict)
mdata.write(snakemake.output[0])
smk_logger.info('Completed writing features')
