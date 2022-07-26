# Make a SLURM Profile
## 1) Make profile from cookiecutter template
```
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
profile_name [slurm]: slurm.my_account
sbatch_defaults []: mail-user=account@email.com mail-type=END,FAIL
cluster_sidecar_help: [Use cluster sidecar. NB! Requires snakemake >= 7.0! Enter to continue...]
Select cluster_sidecar:
1 - yes
2 - no
Choose from 1, 2 [1]:
cluster_name []:
cluster_config_help: [The use of cluster-config is discouraged. Rather, set snakemake CLI options in the profile configuration file (see snakemake documentation on best practices). Enter to continue...]
cluster_config []:
```
## 2) Add conda activation script
Edit the `slurm.my_account/slurm-jobscript.sh` bash script and add `source PATH_TO_CONDA/etc/profile.d/.conda.sh`. Replace `PATH_TO_CONDA` with the path to your conda installation. 

```
#!/bin/bash
# properties = {properties}

source /nethome/$USER/miniconda3/etc/profile.d/conda.sh

{exec_job}
```

## 3) Move the profile to the correct location 
```
mkdir ~/.config/snakemake
mv slurm.my_account/ ~/.config/snakemake/
```

# Running pipeline
```
snakemake --configfile ../config/config.yaml --profile slurm.my_account --use-conda -j 2
```




