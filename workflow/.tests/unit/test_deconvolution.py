import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_deconvolution():
    configuration = ['config_deconwolf_generate_psf','config_deconwolf_use_psf',"config_care",'config_urnet']
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        for config in configuration:
            data_path = PurePosixPath(f".tests/unit/deconvolution/data/{config}")
            expected_path = PurePosixPath(".tests/unit/deconvolution/expected/{config}")

            # Copy data to the temporary workdir.
            shutil.copytree(data_path, workdir)

            # dbg
            print("/gpfs/commons/home/ecordina/ecordina_innovation/Hudson_Test/20240913_BeadsSequencing/deconvolved/deconvolved_Pad2.zarr", file=sys.stderr)

            # Run the test job.
            sp.check_output([
                "python",
                "-m",
                "snakemake", 
                "/gpfs/commons/home/ecordina/ecordina_innovation/Hudson_Test/20240913_BeadsSequencing/deconvolved/deconvolved_Pad2.zarr",
                "-f", 
                "-j1",
                "--keep-target-files",
                "--configfile",
                f"/gpfs/commons/home/ecordina/ecordina_innovation/200/hudson/config/{config}.yaml",
        
                "--use-conda",
                "--directory",
                workdir,
            ])

            # Check the output byte by byte using cmp.
            # To modify this behavior, you can inherit from common.OutputChecker in here
            # and overwrite the method `compare_files(generated_file, expected_file), 
            # also see common.py.
            common.OutputChecker(data_path, expected_path, workdir).check()
