from pyseq import image_analysis as ia
from pre import utils
from os.path import join
from utils import get_cluster




experiment_config = utils.get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])
image_path = join(exp_dir, image_path)

section_name = snakemake.params.section

image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name)

# Start dask cluster
winfo = snakemake.config.resources.dask_worker
p = snakemake.config.resources['partition']
cluster, client = get_cluster(log_dir=None, queue_name = p, **winfo)
print('Client::', cluster.dashboard_link)
ntiles = int(len(im.im.col)/2048)
min_workers = 2*ntiles
max_workers = 2*min_workers
cluster.adapt(minimum = min_workers, maximum=max_workers)
client.wait_for_workers(int(min_workers/2), 60*5)

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
