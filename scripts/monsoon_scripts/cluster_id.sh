#!/bin/bash

#SBATCH --time=30:00
#SBATCH --mem=20G
#SBATCH --job-name=clust

in_file="$1"
out_dir="$2"

# for each fasta file in the input directory:
cd "$in_file"

for file in $(ls *.fasta); do
    # create a directory in out_dir of the name of that file
    mkdir "../$out_dir"/$file

    # echo sbatch commands into the script
    echo "#!/bin/bash" >> "$file".sh
    echo "#SBATCH --time=2:00:00" >> "$file".sh
    echo "#SBATCH --mem=30G" >> "$file".sh
    echo "#SBATCH --job-name=cluster_$file" >> "$file".sh
    echo "module load python/3.latest" >> "$file".sh
    clustering=".././clustering.py -q $file -c kmer --id 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,0.999 -n 200000 -k 9 -o ../$out_dir$file/"
    echo "srun $clustering" >> "$file".sh
    # echo srun and the clustering command into the scrupt

    # run the script
    chmod +x "$file".sh
    sbatch "$file".sh

done


