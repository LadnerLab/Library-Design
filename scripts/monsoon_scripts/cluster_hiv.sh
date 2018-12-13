#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --job-name=do_clust


module load python/3.latest
srun ./clustering.py -q "$1" -o "$2" -c kmer --id 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,0.999 -n 200000 -k 9

