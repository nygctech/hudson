from pyseq import image_analysis as ia
from pre import utils
from os.path import join

experiment_config = utils.get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])
image_path = join(exp_dir, image_path)

section_name = snakemake.params.section

image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name)

# Print out info about section
print('machine::', image.machine)
print('image path::',image_path)
print('section::', section_name)

# Correct Background
print('Correcting background')
print('Pixel group adjustments')
for ch, values in image.config.items(image.machine+'background'):
    print(f'Channel {ch}::',values)

image.correct_background()
image.register_channels2()

#
image.save_zarr(snakemake.params.save_path)

print(image.machine)
