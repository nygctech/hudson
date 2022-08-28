#!/bin/bash

#SBATCH --job-name=mean_intensity
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=6
#SBATCH --mem=36G 
#SBATCH --gres gpu:tesla:1
#SBATCH --time=10:00:00
#SBATCH --output=celltype.log               # Standard output and error log 



snakemake --configfile ../config/config.yaml --profile slurm --use-conda -j 2 --conda-frontend conda

