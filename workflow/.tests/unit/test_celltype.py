import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common
import yaml
import zipfile

def test_celltype():
    
    with TemporaryDirectory(dir='.') as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/celltype/data")
        expected_path = PurePosixPath(".tests/unit/celltype/expected")
        output_path = workdir / 'tables' / 'm1a.csv'

        with zipfile.ZipFile(f"{data_path}/features/m1a.h5mu.zip", 'r') as zip_ref:
            zip_ref.extract('m1a.h5mu', data_path / 'features')
        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Make test config
        test_config_path = common.update_outputdir(workdir)

        # Run the test job.
        with open('test_celltype.out', "w") as outfile:
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
                '--allowed-rules', 'celltype'
            ],
           stderr = outfile)

        exp_files = ['tables/m1a.csv']
        checker = common.OutputChecker(data_path, expected_path, workdir)
        for f in exp_files: 
            assert checker.compare_files(expected_path / f, workdir / f)
