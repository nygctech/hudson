import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common
import yaml


def test_picasso():

    with TemporaryDirectory(dir = '.') as tmpdir:

        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/picasso/data")
        output_path = workdir / 'final_zarr' / 'picasso_m1a.yaml'

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Make test config
        test_config_path = common.update_outputdir(workdir)
        

        # Run the test job.
        with open('test_picasso.out', "w") as outfile:
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
                '--allowed-rules', 'picasso'
            ], 
            stdout = sp.PIPE, stderr = outfile)


        with open(output_path) as f:
            picasso_params = yaml.safe_load(f)
            
            # Test structure of picasso parameter yaml file
            # For each cycle a list of channels and sinks should exist, the list may be empty
            # If a sink is listed a channel must be listed too, and vice versa
            # If a sink/channel is listed background and alpha parameters must be listed
            # The background/alpha parameter matrix will be nchannel x nsinks in size
    
            for cycle, params in picasso_params.items():
                chs = params.get('channels', None)
                sinks = params.get('sinks', None)
                assert chs is not None
                assert sinks is not None
                nsinks = len(sinks)
                nparams = len(chs)
                if nsinks > 0:
                    for s in sinks:
                        assert s in chs
                        assert isinstance(s, int)
                    bg = params.get('background', None)
                    a = params.get('alpha', None)
                    assert bg is not None
                    assert a is not None
                    assert len(bg) == nparams
                    assert len(a) == nparams
                    for i in range(nparams):
                        assert len(bg[i]) == nsinks
                        assert len(a[i]) == nsinks
