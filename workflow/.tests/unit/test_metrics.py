import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common
import yaml


def test_metrics():
    
    with TemporaryDirectory(dir='.') as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/metrics/data")
        expected_path = PurePosixPath(".tests/unit/metrics/expected")
        output_path = workdir / 'celltype' / 'm1a.h5mu'
        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Make test config
        test_config_path = common.update_outputdir(workdir)


        with open('test_metrics.out', "w") as outfile:
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
                '--allowed-rules', 'segmentation', 'feature_extraction', 'celltype'
            ],
           stderr = outfile)

        exp_files = ['features/m1a.csv', 'summary_m1a.yaml']
        shutil.copy(workdir/ 'summary_m1a.yaml', '/commons/groups/innovation/sarah/hudson/workflow/summary_m1a.yaml')
        checker = common.OutputChecker(data_path, expected_path, workdir)
    
        for f in exp_files: 
            assert checker.compare_files(expected_path / f, workdir / f)