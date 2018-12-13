#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=1500G
#SBATCH -c 1
#SBATCH -o hiv_validate.txt_coll

module load python/2.7.5
echo "$1"
echo "$2"
srun ./validate_design.py -g "$1" -r "$2" -k 7,8,9,10,11,12,13  --makemap -o "$1"_map
