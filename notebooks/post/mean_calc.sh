#!/bin/bash
#SBATCH --job-name=mean_calc                  # Job name
#SBATCH --partition=pe2                     # Partition Name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsingh@nygenome.org     # Where to send mail
#SBATCH --mem=16G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=01:00:00                       # Time limit
#SBATCH --output=mean_calc.log               # Standard output and error log 

#conda activate spatial
python -u test_run.py mean_calc

