#!/bin/bash
#SBATCH --mail-type=end
#SBATCH --mail-user=seandavi@gmail.com
#SBATCH --partition=ccr
#SBATCH --time=72:00:00
module load python/3.4.3
snakemake -j 400 --drmaa -k --ri --stats snakemake.stats
