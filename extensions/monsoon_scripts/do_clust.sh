#!/bin/bash
#SBATCH --time=2-23:05:10
#SBATCH --mem=1400G
#SBATCH -o clust_redone.txt
#SBATCH --open-mode=append
#SBATCH --job-name=clust_900000
#SBATCH --array=0-2
#SBATCH -c 15
cores=15

if [ "$#" -lt 2 ]; then
    echo Not enough arguments.
    echo Usage:
    echo sbatch do_clust.sh input_dir logfile_name

    exit 1
fi

input_dir="$1"
tracking_log="$3"


cd "$input_dir"

files=($(ls *fasta* | grep -v "out"))
file=${files[$SLURM_ARRAY_TASK_ID]}

echo "$file"'|'"$SLURM_ARRAY_JOB_ID"_"$SLURM_ARRAY_TASK_ID" >> "$tracking_log"

if [ ! -e "$file"_out* ]; then
    echo Doing "$file"
    srun ../../kmer_oligo -q $file -i 1 -o "$file"_out -x 9 -y 24 -t "$cores"
fi
