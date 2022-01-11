#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=2gb
#SBATCH --time=120:00:00
#SBATCH --account=cschmal1
#SBATCH --mail-user=cschmal1@uni-koeln.de
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --output=/home/cschmal1/phQTL/output/Perlstein-%j
#SBATCH --error=/home/cschmal1/phQTL/output/Perlstein-%j
#SBATCH --array=1-307

module load R
export R_LIBS_USER=$HOME/R/3.1.1
Rscript --vanilla --verbose /home/cschmal1/phQTL/scripts/Perlstein_QTLmapping_cheops.R ${SLURM_ARRAY_TASK_ID}


