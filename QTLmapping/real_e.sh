#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=12gb
#SBATCH --time=120:00:00
#SBATCH --account=jgrossb1
#SBATCH --output=/scratch/jgrossb1/phosphoQTL/eQTL/out/o-%j
#SBATCH --error=/scratch/jgrossb1/phosphoQTL/eQTL/out/e-%j
#SBATCH --mail-user=jgrossb1@uni-koeln.de
#SBATCH --mail-type=FAIL

Rscript --vanilla --verbose real_e.R

