from picasso.nn_picasso import PICASSOnn
from pre import image_analysis as ia

image = ia.get_HiSeqImages(snakemake.input[0])

stains_config = snakemake.config.get('stains')
sink_source_ch = snakemake.config.get('unmixing', {610:[558], 740:[687]})

channels = list(image.im.channels)
cycles = list(image.im.cycles)

cy_mm = {}

# See if background or autofluorescence image exists
AF = ()
for cy, ch_stains in stains_config.items():
    for ch, marker in ch_stains.items():
        if marker in ['background', 'autofluorescence', 'af', 'AF', 'bg', 'BG']:
            AF = (cy,ch)

# Find which channels to unmix
for cy, ch_stains in stains_config.items():
    si_so = {}                                                                  # sink channel : [source channels]
    all_ch = []                                                                 # all channels that need to be used for unmixing

    # Get dictionary of sink channels that need source channels removed
    for ch in ch_stains.keys():                                                 # loop over channels used in cycle
        for source_ch in sink_source_ch.get(ch,[]):                             # loop over source channels if a channel has spillover fluorescence
            if source_ch in ch_stains.keys():                                   # add source channel to sink_source dict (si_so)
                si_so.setdefault(sink_ch),[]).append(source_ch)
                all_ch.append(all_ch)
        all_ch =+ si_so.keys()
        all_ch = np.unique(all_ch)

    # Making mixing matrix for cycle
    mm = np.array(len(all_ch), len(si_so))
    for c, si, so_ in enumerate(si_so.items()):
        si_ind = np.where(all_ch == si)
        mm[si_ind, c] = 1
        for so in so_:
            so_ind = np.where(all_ch == so)
            mm[so_ind, c] = -1
    cy_mm[cy] = [mm, list(all_ch), list(si_so.keys())]                          # dict with cycle keys and value = list with mixing matrix, all images, and all sink image


# Unmix images and stack by stain instead of cycle > channel
stain_stack = []
stain_name = []
for cy , ch_stains in stains_config.items():
    print(f'Unmixing cycle {cy}')

    mm, all_ch, sinks = cy_mm.get(cy,[None, None, []])
    if len(sinks) > 0:
        print('Removing fluorescence from channels', *sinks)
        model = PICASSOnn(cy_mm[cy][0])
        mixing_params = model.train_loop(images = images.sel(cycle=cy, channel=all_ch)
        print('Mixing params:', mixing_params)
        unmixed_ims = model.unmix_images(images = images.sel(cycle=cy, channel=all_ch)
        unmixed_ims = xr.DataArray(unmixed_ims, dim=['channel','row','col'], coords = {'channel':sinks})

    for ch, marker = ch_stains.items():
        if ch in sinks:
            stain_stack.append(unmixed_ims.sel(channel=ch))                     # add unmixed sink
        else:
            stain_stack.append(image.im.sel(cycle=cy, channel=ch))              # add channel that did not have any spillover
        stain_name.append(marker)

# Save image
unmixed = xr.concat(stain_stack, dim = 'marker', coordinates = {'marker':stain_name})
delayed_store = unmixed.to_dataset().to_zarr(snakemake.params.save_path, compute = False)
