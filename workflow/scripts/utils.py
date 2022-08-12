
def get_cluster(manager, **winfo):

    import re

    assert manager.upper() in ['SLURM']

    if re.search(manager, 'SLURM', re.IGNORECASE):
        from dask_jobqueue import PBSCluster
        cluster = SLURMCluster(**winfo)

    return cluster
