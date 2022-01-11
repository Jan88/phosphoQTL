#!/bin/bash -l
#SBATCH -J=pQTL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=120:00:00
#SBATCH --account=UniKoeln
#SBATCH --mail-user=jgrossb1@uni-koeln.de
#SBATCH --array=1-133
#SBATCH --share

module load R/3.3.3_*
R --vanilla -f permutations_pt.R --args ${SLURM_ARRAY_TASK_ID}



