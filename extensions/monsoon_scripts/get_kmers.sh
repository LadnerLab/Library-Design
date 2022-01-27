#!/bin/bash
#SBATCH --mem=3M
#SBATCH --array=0-102
#SBATCH --time=00:30
#SBATCH -o kmer_count_orig_clusters.txt
#SBATCH --open-mode=append

module load python/3.latest

files=$(ls "$1"/*.fasta)
file=${files[$SLURM_ARRAY_TASK_ID]}
echo $file
#srun count_kmers.py $file
