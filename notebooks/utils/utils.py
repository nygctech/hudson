import pathlib
import configparser
from os.path import join
from os import getcwd, makedirs
from dask_jobqueue import SLURMCluster

#dask.config.set({'temporary-directory': '/scratch'})
def get_cluster(engine = 'SLURM', **kwargs):
    """ Make dask cluster w/ workers = 2 cores, 32 G mem, and 1 hr wall time.

        Parameters:
            engine: name of workload manager. Only SLURM currently supported.


        Returns:
            cluster: dask cluster
    """

    assert engine in ['SLURM'], f'{engine} is not currently supported as a workload manager'

    if kwargs.get('log_directory', None) is None:
        kwargs['log_directory'] = join(getcwd(),'dask_logs')
        makedirs(kwargs['log_directory'], exist_ok=True)

    if engine == 'SLURM':
        cluster = SLURMCluster(**kwargs)
                # queue = queue_name,
                # cores = 2,
                # memory = '32G',
                # walltime='1:00:00',
                # log_directory=log_dir,
                # extra=["--lifetime", "55m", "--lifetime-stagger", "4m"])

    print('Cluster dashboard link::', cluster.dashboard_link)

    return cluster
