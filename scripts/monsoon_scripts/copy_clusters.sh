#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --job-name=copy_clusters

module load python/3.latest

input_dir="$1"
cluster_min="$2"
cluster_max="$3"

if [[ $# -lt 3 ]]; then
    echo "Usage:"
    echo ./copy_clusters input_dir cluster_min cluster_max
    exit 1
fi
cluster_items=($(ls "$input_dir" | grep -v "out"))
mkdir "$1"/clusters_"$cluster_max"
for item in "${cluster_items[@]}"; do
    item="$input_dir"/"$item"
    count=$(./count_kmers.py "$item" | cut -d '|' -f 2)
    if [ "$count" -ge "$cluster_min" ] && [ "$count" -lt "$cluster_max" ]; then 

        new_name=$( ls "$item" | cut -d '/' -f 2)
        cp "$item" "$input_dir"/clusters_"$cluster_max"/"$new_name"
    fi
done  
