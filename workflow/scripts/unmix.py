from picasso.nn_picasso import PICASSOnn
import xarray as xr
import numpy as np
from pathlib import Path
import dask

# Open image from zarr store
image_path = Path(snakemake.input[0])
im_name = image_path.stem
image = xr.open_zarr(image_path).to_array()
image = image.squeeze().drop_vars('variable').rename(im_name)

# Get config data
stains_config = snakemake.config.get('stains')
sink_source_ch = snakemake.config.get('unmixing', {610:[558], 740:[687]})

# See if background or autofluorescence image exists
AF = ()
for cy, ch_stains in stains_config.items():
    for ch, marker in ch_stains.items():
        if marker in ['background', 'autofluorescence', 'af', 'AF', 'bg', 'BG']:
            AF = (cy,ch)


# Find which channels to unmix
cy_mm = {}
for cy, ch_stains in stains_config.items():
    si_so = {}                                                                  # sink channel : [source channels]
    all_ch = []                                                                 # all channels that need to be used for unmixing

    # Get dictionary of sink channels that need source channels removed
    for ch in ch_stains.keys():                                                 # loop over channels used in cycle
        for source_ch in sink_source_ch.get(ch,[]):                             # loop over source channels if a channel has spillover fluorescence
            if source_ch in ch_stains.keys():                                   # add source channel to sink_source dict (si_so)
                si_so.setdefault(ch,[]).append(source_ch)
                all_ch.append(source_ch)
    all_ch =+ si_so.keys()
    all_ch = np.unique(all_ch)

    # Making mixing matrix for cycle
    mm = np.zeros((len(all_ch), len(si_so)))
    for c, si, so_ in enumerate(si_so.items()):s
        si_ind = np.where(all_ch == si)
        mm[si_ind, c] = 1
        for so in so_:
            so_ind = np.where(all_ch == so)
            mm[so_ind, c] = -1

    cy_mm[cy] = [mm, list(all_ch), list(si_so.keys())]                        # dict with cycle keys and value = list with mixing matrix, all images, and all sink image


# Unmix images and stack by stain instead of cycle > channel
stain_stack = []
stain_name = []
for cy , ch_stains in stains_config.items():
    mm, all_ch, sinks = cy_mm.get(cy)

    if len(sinks) > 0:
        print(f'Unmixing cycle {cy}')
        print('Removing fluorescence from channels', *sinks)
        model = PICASSOnn(cy_mm[cy][0])

        # Fit alpha and background
        for i in model.train_loop(images = image.sel(cycle=cy, channel=all_ch)):
            pass
        print('Mixing params:')
        coords = {'channel':all_ch, 'sink':sinks}
        alpha = xr.DataArray(model.mixing_parameters[0,:,:], name = 'alpha' dims = dims, coords = coords)
        bg = xr.DataArray(model.mixing_parameters[1,:,:], name = 'background', dims = dims, coords = coords)
        print(alpha)
        print(bg)

        # Unmix Images
        ims = image.sel(cycle=cy, channel=all_ch).stack(dims = ['row','col']).T
        ch_stack = []
        for sink_ch in sinks:
            unmixed_im = ((ims - bg.sel(sink=sink_ch)) @ alpha.sel(sink=sink_ch)).unstack().drop_vars('sink')
            ch_stack.append(unmixed_im)
        unmixed_ims = xr.concat(ch_stack, dim='channel').assign_coords({'channel':sinks})

    # Save as markers instead of cycle & channel
    for ch, marker = ch_stains.items():
        if ch in sinks:
            stain_stack.append(unmixed_ims.sel(channel=ch))                     # add unmixed sink
        else:
            stain_stack.append(image.im.sel(cycle=cy, channel=ch))              # add channel that did not have any spillover
        stain_name.append(marker)

# Save image
unmixed = xr.concat(stain_stack, dim='marker').assign_coords({'marker':stain_name})
delayed_store = unmixed.to_dataset().to_zarr(snakemake.params.save_path, compute = False)
dask.compute(delayed_store)
