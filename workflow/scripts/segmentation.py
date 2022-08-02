import xarray as xr
import cellpose
from pathlib import Path
import io
import imageio


#Open image from zarr store
image_path = Path(snakemake.input[0])
im_name = image_path.stem
image = xr.open_zarr(image_path).to_array()
image = image.squeeze().drop_vars('variable').rename(im_name)

# segment
logger = io.logger_setup()
use_GPU = int(snakemake.params.get('use_GPU'))
model = models.CellposeModel(gpu=use_GPU,model_type='TN2')
#model = models.CellposeModel(model_type='TN2')
# Remove once priors steps in pipe
#one_z_plane = image.sel(obj_step = 8498, channel = 558, cycle=1)
arr = image.sel(snakemake.config.get('segmentation'))
channels = [0,0]
masks, flows, styles = model.eval(arr.values, diameter=None, channels=channels, cellprob_threshold= -6, flow_threshold= 1000)
imageio.imwrite(snakemake.output[0],masks)

