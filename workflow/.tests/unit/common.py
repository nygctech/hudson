"""
Common code for unit testing of rules generated with Snakemake 7.12.0.
"""

from pathlib import Path
import subprocess as sp
import os
import yaml
# import string
# import random



# # Copied from stack overflow
# def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
#     return ''.join(random.choice(chars) for _ in range(size))

def update_outputdir(workdir):
    # TODO update with path from remote data repo
    config_path = Path(os.path.abspath(__file__)).parent / 'pytest_config.yaml'

    with open(config_path) as f:
        config = yaml.safe_load(f)

    config.update({'output_directory': str(workdir)})
    
    test_config_path = workdir / 'pytest_config.yaml'
    with open(test_config_path, 'w') as f:
        f.write(yaml.dump(config))
    
    return test_config_path

class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if str(f).startswith(".snakemake"):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):
        # KP Modified to return true if files are the same and catch errors if not
        try: 
            sp.check_output(["cmp", generated_file, expected_file])
            return True
        except Exception as error:
            print(type(error).__name__, "–", error)
            print(error.output)
            return False
 
