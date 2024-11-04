import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common
import yaml
import xarray as xr
import numpy as np 

def test_saverae():
    

    with TemporaryDirectory(dir = '.') as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/saveraw/data")
        expected_path = PurePosixPath(".tests/unit/saveraw/expected")
        output_path = workdir / 'raw_zarr' / 'm1a.zarr'

        shutil.copytree(data_path, workdir)
        # Make test config
        test_config_path = common.update_outputdir(workdir)
        # Run the test job.
        with open('test_saveraw.out', "w") as outfile:
            sp.check_output([
                "python",
                "-m",
                "snakemake", 
                output_path,
                "-f", 
                "-j1",
                "--keep-target-files",
                "--configfile",
                test_config_path,
                "--use-conda",
                '--allowed-rules', 'saveraw'
            ],
            stderr = outfile)

        exp_files = ['raw_zarr/m1a.zarr']
        with open(test_config_path) as f:
            exp_config = yaml.safe_load(f)

        cycles = exp_config['markers']
        data = xr.open_zarr(workdir / exp_files[0])

        for cycle in list(cycles.keys()):
            assert cycle in data.cycle.values
        
        markers = []
        for key, value in cycles.items():
            if isinstance(value, dict):  # Check if the value is a dictionary
                markers.extend(value.keys()) 
        
        for marker in np.unique(markers):
            assert marker in data.channel.values
