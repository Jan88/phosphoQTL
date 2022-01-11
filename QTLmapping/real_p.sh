#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8 
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --time=120:00:00
#SBATCH --account=jgrossb1
#SBATCH --output=/scratch/jgrossb1/phosphoQTL/pQTL/out/o-%j
#SBATCH --error=/scratch/jgrossb1/phosphoQTL/pQTL/out/e-%j
#SBATCH --mail-user=jgrossb1@uni-koeln.de
#SBATCH --mail-type=FAIL

module load R/3.3.3_*
R --vanilla -f real_p.R
