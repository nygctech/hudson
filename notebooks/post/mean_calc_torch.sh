#!/bin/bash

#SBATCH --job-name=mean_intensity
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=6
#SBATCH --mem=36G 
#SBATCH --gres gpu:tesla:1
#SBATCH --time=5:00:00
#SBATCH --output=mean_calc_gpu.log               # Standard output and error log 

python cellpose.py 


