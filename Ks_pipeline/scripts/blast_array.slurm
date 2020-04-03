#!/bin/bash
#SBATCH --mail-user=jacob.steenwyk@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=20G
#SBATCH --array=0-9
#SBATCH --output=slurm_blast_array_%A_%a.out

bash x0${SLURM_ARRAY_TASK_ID}
