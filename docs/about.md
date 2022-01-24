---
layout: default
title: About
permalink: /about/
---
Current Version: v1.3.0

### Requirements
- Slurm-Based HPC Cluster
- Python 3.4 or greater


### Input
A set of clusters (1 or more), each containing protein sequences. Clusters can be generated through any means,
including MUSCLE, [uniprotfasta_to_clusters.py](scripts/uniprotfasta_to_clusters.py), etc.

### Output
A set of k-mer oligonucleotides that can be synthesized for use in serological detection of viral exposure history.

### Step 1: Grouping Clusters of Similar Kmer-Sizes
**Note** This step can be skipped if desired, you will however make your cluster administrator
happier if you group clusters, as resource demands of large clusters ( > 200,000 k-mer epitopes ) are far greater than the
demands of smaller clusters (0-50,000 k-mer epitopes).

Also note that cluster size bounds take the form [ a, b ), a < b.
Experimentally, the following cluster groupings have been successful in
testing:
* [ 1, 50,000 )
* [ 50,000, 100,000 )
* [ 100,000, 200,000 )
* [ 200,000, 1,000,000 )


For each of the above thresholds, run the following command from the directory outside
of your clusters:
```
./copy_clusters.sh a b
```
where a is the lower bound (inclusive), and b is the upper bound (exclusive).
For each of the above thresholds, the original cluster directory now contains subdirectories
of the form clusters_a_b.

### Step 2: Create the Oligos
**Note** for information concerning to the usage of this script, please see the README found at
[this link.](https://github.com/LadnerLab/C-KmerOligo)

Edit the `srun` line in [do_clust.sh](scripts/monsoon_scripts/do_clust.sh) with the parameters
you prefer.

For each of the directories containing groups of clusters, do the following:
- Obtain the number of clusters are contained in the grouped directory:

```
ls base_dir/clusters_a_b/*.fasta | wc -l
```
- Update the #SBATCH lines in [do_clust.sh](scripts/monsoon_scripts/do_clust.sh) with the
  parameters appropriate for the size of the clusters. The array line should be updated to
  reflect the number of fasta files in the directory.

- Run the script:
```
sbatch do_clust.sh base_dir/clusters_a_b logfile_name.txt
```


### Step 3: Combine output files and remove duplicated Oligos
- Combine the output files:
```
awk 1 base_dir/clusters_*/*_out_R_1 >> combined_file.fasta
```

Often times there will be a few oligos that are not unique in this combined file,
remove them with:
```
remove_duplicate_oligos.py combined_file.fasta
```
This will produce the file combined_file.fasta_unique

### Step 4 (Optional): Produce map file
A map file contains entries that pair every produced oligonucleotide
with every sequence from the input containing a shared 7-mer.

```
sbatch validate_maponly.sh combined_file.fasta_unique unclustered_complete_input.fasta
```
