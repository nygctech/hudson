## Make a profile
# Make profile from cookiecutter template
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
# Move the profile to the correct location 
```
mkdir ~/.config/snakemake
mv slurm.my_account/ ~/.config/snakemake/
```

# Running pipeline
```
conda activate pyseq_pipe
snakemake --configfile ../config/config.yaml --profile slurm.my_account -j 2
```




