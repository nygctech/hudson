
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/nethome/kpandit/miniconda3/envs/pyseq_pipe/lib/python3.9/site-packages', '/gpfs/commons/home/kpandit/hudson/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95r\x06\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c5/nethome/kpandit/pyseq-image/src/demo/logs/config.cfg\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x10\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x16)}\x94\x8c\x05_name\x94h\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cI/nethome/kpandit/hudson/results/20210323_4i4color/processed_zarr/m1a.zarr\x94a}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8c@/nethome/kpandit/hudson/results/20210323_4i4color/processed_zarr\x94\x8c\x03m1a\x94e}\x94(h\x0c}\x94(\x8c\tsave_path\x94K\x00N\x86\x94\x8c\x07section\x94K\x01N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bh8h4h:h5ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c1/nethome/kpandit/hudson/results/20210323_4i4color\x94h5e}\x94(h\x0c}\x94(\x8c\noutput_dir\x94K\x00N\x86\x94\x8c\x08sections\x94K\x01N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94b\x8c\noutput_dir\x94hI\x8c\x08sections\x94h5ub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01M\xe8\x03M\xe8\x03\x8c\x08/scratch\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07disk_mb\x94K\x03N\x86\x94\x8c\x06tmpdir\x94K\x04N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bhcK\x01heK\x01hgM\xe8\x03hiM\xe8\x03hkh`ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x14experiment_directory\x94\x8c&/nethome/kpandit/pyseq-image/src/demo/\x94\x8c\x10output_directory\x94\x8c /nethome/kpandit/hudson/results/\x94\x8c\tresources\x94}\x94(\x8c\tpartition\x94\x8c\x03pi3\x94\x8c\x04time\x94\x8c\x071:00:00\x94\x8c\x03mem\x94\x8c\x0316G\x94\x8c\x0bdask_worker\x94}\x94(\x8c\x05cores\x94K\x02\x8c\x05extra\x94]\x94(\x8c\n--lifetime\x94\x8c\x0355m\x94\x8c\x12--lifetime-stagger\x94\x8c\x024m\x94e\x8c\x06memory\x94\x8c\x0332G\x94\x8c\x08walltime\x94\x8c\x071:00:00\x94uuu\x8c\x04rule\x94\x8c\x0cfix_lighting\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c2/gpfs/commons/home/kpandit/hudson/workflow/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/gpfs/commons/home/kpandit/hudson/workflow/scripts/preprocess.py';
######## snakemake preamble end #########
#from pyseq import image_analysis as ia
from pre import utils
from pre import image_analysis as ia
from os.path import join
from utils.utils import get_cluster




experiment_config = utils.get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
image_path = snakemake.config.get('image_path',experiment_config['experiment']['image path'])
image_path = join(exp_dir, image_path)

section_name = snakemake.params.section

image = ia.get_HiSeqImages(image_path = image_path, common_name = section_name)

# Start dask cluster
default_winfo = {'cores':2, 'memory':'32G', 'walltime':'1:00:00'}
winfo = snakemake.config.get('resources',{}).get('dask_worker', default_winfo)
p = snakemake.config['resources']['partition']
cluster = get_cluster(log_dir=None, queue_name = p, **winfo)
ntiles = int(len(image.im.col)/2048)
min_workers = max(1,2*ntiles)
max_workers = 2*min_workers

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
