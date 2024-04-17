import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common
import yaml


def test_feature_extraction():
    
    with TemporaryDirectory(dir='.') as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        #workdir = Path('./20240417/workdir')
        data_path = PurePosixPath(".tests/unit/feature_extraction/data")
        expected_path = PurePosixPath(".tests/unit/feature_extraction/expected")
        output_path = workdir / 'features' / 'm1a.h5mu'

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Make test config
        test_config_path = common.update_outputdir(workdir)
        print(f'{test_config_path=}')

        # Increase SLURM resources for test
        with open(test_config_path) as f:
             config = yaml.safe_load(f)
        config.update({'resources': 
                        {'dask_worker':
                          {'cores': 1,
                           'memory': '8G',
                           'manager': 'SLURM',
                           'log_directory': 'dask_logs'
                      }}})
        with open(test_config_path, 'w') as f:
             f.write(yaml.dump(config))

        # dbg
        print(f'{workdir=}')
        print(f'{data_path=}')
        print(f'{expected_path=}')
        print(f'{output_path=}')


        # Run the test job.
        with open('test_feature_extraction.out', "w") as outfile:
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
                #"-n",
                # "--directory",
                # workdir,
            ],
           stderr = outfile)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        #common.OutputChecker(data_path, expected_path, workdir).check()

        #common.OutputChecker(data_path, expected_path, workdir).check()
        exp_files = ['features/m1a.h5mu',
                    ]
        checker = common.OutputChecker(data_path, expected_path, workdir)
    	
        for f in exp_files: 
            assert checker.compare_files(expected_path / f, workdir / f)
