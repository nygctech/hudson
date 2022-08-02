# import
import xarray as xr
from pathlib import Path

class snakemake()
    input = ['path to input']
    output = ['path to output']


# Open image from zarr store
image_path = Path(snakemake.input[0])
im_name = image_path.stem
image = xr.open_zarr(image_path).to_array()
image = image.squeeze().drop_vars('variable').rename(im_name)

# segmentation


# save it to
imageio.imwrite(mask, snakemake.output[0])
