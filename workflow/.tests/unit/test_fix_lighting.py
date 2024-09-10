import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_saveraw():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"

        data_path = PurePosixPath(".tests/unit/fix_lighting/data")
        
        expected_path = PurePosixPath(".tests/unit/fix_lighting/expected")

        output_path = workdir/ 'processed_zarr' / 'm1a.zarr'

        shutil.copytree(data_path, workdir)
        
        test_config_path = common.update_outputdir(workdir)


        with open('test_fix_lighting.out', "w") as outfile:
            sp.run([
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
         #       "--unlock",
         #       "--directory",
         #       workdir,
            ], 
            stdout = sp.PIPE, stderr = outfile)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        exp_files = ['processed_zarr/m1a.zarr',
                    ]
        checker = common.OutputChecker(data_path, expected_path, workdir)
    	
        for f in exp_files: 
            assert checker.compare_files(expected_path / f, workdir / f)



