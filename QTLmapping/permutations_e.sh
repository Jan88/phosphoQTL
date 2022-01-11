#!/bin/bash -l
#SBATCH -J=eQTL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=2GB
#SBATCH --time=120:00:00
#SBATCH --account=UniKoeln
#SBATCH --mail-user=jgrossb1@uni-koeln.de
#SBATCH --array=1-543
#SBATCH --share

Rscript --vanilla --verbose permutations_e.R ${SLURM_ARRAY_TASK_ID}


