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
        # tmpdir = Path('./20240416')
        # workdir = tmpdir / "workdir"
        data_path = PurePosixPath(".tests/unit/unmix/data")
        expected_path = PurePosixPath(".tests/unit/unmix/expected")
        output_path = workdir / 'final_zarr' / 'm1a.zarr'

        print(f'{workdir=}')
        print(f'{data_path=}')
        print(f'{expected_path=}')
        print(f'{output_path=}')


        # Copy data to the temporary workdir.
        # shutil.copytree(data_path, workdir)

        # Make test config
        test_config_path = common.update_outputdir(workdir)
        print(f'{test_config_path=}')

        # dbg
        # print("/gpfs/commons/groups/innovation/sarah/Hudson_Test/ome_zarr_test_2/20210323_4i4color/final_zarr/m1a.zarr", file=sys.stderr)

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
                "--profile",
                "slurm",
                #'-n',
                # "--directory",
                # workdir,
                # "--unlock",
                # "--rerun-incomplete",
            ],
            #stdout = sp.PIPE, 
            stderr = outfile)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        #common.OutputChecker(data_path, expected_path, workdir).check()
        exp_files = ['final_zarr/m1a.zarr/.zattrs',
                     'final_zarr/m1a.zarr/0/.zarray'
                    ]
        checker = common.OutputChecker(data_path, expected_path, workdir)
        for f in exp_files: 
            assert checker.compare_files(expected_path / f, workdir / f)

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
            
