#!/bin/bash
#SBATCH --mail-user=jacob.steenwyk@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=20G
#SBATCH --output=01_slurm.out

SHORT1=$(cat file.ctl | awk '{print $1}')
CDS1=$(cat file.ctl | awk '{print $2}')
PEP1=$(cat file.ctl | awk '{print $3}')

SHORT2=$(cat file.ctl | awk '{print $4}')
CDS2=$(cat file.ctl | awk '{print $5}')
PEP2=$(cat file.ctl | awk '{print $6}')

bash ./scripts/determine_Ks_for_every_gene_between_two_species_pt1.bash $SHORT1 $CDS1 $PEP1 $SHORT2 $CDS2 $PEP2


