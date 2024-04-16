import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_fix_lighting():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/fix_lighting/data")
        expected_path = PurePosixPath(".tests/unit/fix_lighting/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("/gpfs/commons/groups/innovation/sarah/Hudson_Test/ome_zarr_test_2/20210323_4i4color/processed_zarr/m1a.zarr", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "/gpfs/commons/groups/innovation/sarah/Hudson_Test/ome_zarr_test_2/20210323_4i4color/processed_zarr/m1a.zarr",
            "-f", 
            "-j1",
            "--keep-target-files",
            "--configfile",
           "/gpfs/commons/groups/innovation/sarah/Hudson_Test/hudson_omezarr/config/config.yaml",
    
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
