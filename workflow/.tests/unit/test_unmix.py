import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common
import yaml

def test_unmix():
    

    with TemporaryDirectory(dir = '.') as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/unmix/data")
        expected_path = PurePosixPath(".tests/unit/unmix/expected")
        output_path = workdir / 'final_zarr' / 'm1a.zarr'
        picasso=workdir/'final_zarr'/"picasso_m1a.yaml"
        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        os.system(f"touch {picasso}")
        # Make test config
        test_config_path = common.update_outputdir(workdir)
        
        # Run the test job.
        with open('test_unmix.out', "w") as outfile:
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
                '--allowed-rules', 'unmix'
            ],
            stderr = outfile)

        exp_files = ['final_zarr/m1a.zarr/.zattrs',
                     'final_zarr/m1a.zarr/0/.zarray'
                    ]
        # Check channel, cycle, marker name matches
        with open(test_config_path) as f:
            exp_config = yaml.safe_load(f)
        with open(workdir / exp_files[0]) as f:
            zattrs = yaml.safe_load(f)            
        out_markers = zattrs['omero']['images'][0]['pixels']['channels']
        for m in out_markers:
            out_name = m['name']
            out_ch = int(m['emission_wavelength'])
            out_cy = int(m['fluor'].split(' ')[1])
            assert out_name == exp_config['markers'][out_cy][out_ch]
