#from pyseq import image_analysis as ia

import os
import sys

print('Python PATH')
os.system('which python')
print('executing python', sys.executable)



from pre import utils
from pre import image_analysis as ia
from os.path import join
from utils.utils import get_cluster
from dask.distributed import Client

experiment_config = utils.get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])
image_path = join(exp_dir, image_path)

section_name = snakemake.params.section

image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name)

# Start dask cluster
# specify default worker options in ~/.config/dask/jobqueue.yaml
winfo = snakemake.config.get('resources',{}).get('dask_worker',{})
cluster = get_cluster(**winfo)
print(cluster.dashboard_link)
print(cluster.job_script())
ntiles = int(len(image.im.col)/2048)
min_workers = max(1,2*ntiles)
max_workers = 2*min_workers
print(min_workers, max_workers)

# Print out info about section
print('machine::', image.machine)
print('image path::',image_path)
print('section::', section_name)

# Start computation
with Client(cluster) as client:

	cluster.adapt(minimum = min_workers, maximum=max_workers)
	client.wait_for_workers(int(min_workers/2), 60*5)


	# Correct Background
	print('Correcting background')
	print('Pixel group adjustments')
	for ch, values in image.config.items(image.machine+'background'):
    		print(f'Channel {ch}::',values)

	image.correct_background()
	image.register_channels2()

#
	image.save_zarr(snakemake.params.save_path)
