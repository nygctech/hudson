
def get_cluster(manager='SLURM', **winfo):

    import re

    assert manager.upper() in ['SLURM']

    if re.search(manager, 'SLURM', re.IGNORECASE):
        from dask_jobqueue import SLURMCluster
        cluster = SLURMCluster(**winfo)

    return cluster
