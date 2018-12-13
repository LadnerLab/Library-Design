#!/bin/bash
#SBATCH --time=01:35:00
#SBATCH --mem=75G
#SBATCH -c 4 
#SBATCH --job-name=kmer_single

srun kmer_oligo -q "$1" -x 9 -y 24 -t 4 -o "$1"_out
