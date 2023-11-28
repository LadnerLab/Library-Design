---
layout: default
title: Tutorial
permalink: /tutorial/
---
# Tutorial

## Overview

This tutorial will walk through the steps to design a library of unique peptides for all of the proteins assigned to the Poxviridae family downloaded from UniProt.

## Inputs

The only required input is a fasta-formatted file containing a set of proteins.

## Outputs

There is one required output, a fasta-formatted file containing the unique, representative sequences from the input. 

There is one optional output, a tab-delimited map file, which relates all removed sequences to the representatives to which they are identical (one line per unique representative).

## Installation

- Visit [this link](https://github.com/LadnerLab/ProteinOligoLibrary/blob/master/protein_oligo_library.py) to download the Protein Oligo Library file and move it to the directory where the scripts you intend to run are located, or place the file in your PYTHONPATH. For example, if you have Anaconda installed, you can move the library to your base environment by putting it here: "~/conda_directory_name/lib/python3.x/site-packages/".
- The installation for the scripts used in this tutorial can be found in each of their respective tabs

### Tutorial/Use

Generate a fasta-formatted file containing your target proteins of interest. For this example, we pulled all of the protein sequences (from Uniprot) for the Poxviridae family.
Remove sequences shorter than the desired output peptide length [subset_fasta_by_length.py]. In this tutorial we will be designing 30AA long peptides.

```
subset_fasta_by_length.py -f poxviridae_unaligned.fasta  -o poxviridae_unaligned_30AA.fasta -s 30
```

(Optional) Remove identical protein sequences, including those that are a subset of another target [onehundredreps.py,one_hundred_reps].

Command (Python version):
```
onehundredreps.py -f poxviridae_unaligned_30AA.fasta -r poxviridae_unaligned_30AA_100rep.fasta -m 100_rep_map.txt
```
Command (C version, using 2 threads):
```
one_hundred_reps \
  -f poxviridae_unaligned_min30AA.fasta \
  -n 2 \
  -m 100_rep_map.txt
```

Make clusters folder to output resulting cluster files into

Command (Linux):
```
mkdir clusters
```

Generate clusters of similar sequences [UCLUST, CD-HIT]. We generally recommend targeting a cluster similarity of 65%-75%.

Command:
```
usearch -cluster_fast poxviridae_unaligned_30AA_100rep.fasta \
    -id 0.70 \
    -sort length \
    -clusters ./clusters/id_70_
```

Make SW\_SC folder within clusters and change working directory to SW\_SC folder

Command:
```
mkdir clusters/SW_SC
cd clusters/SW_SC
```

For each cluster, use a combined sliding window/set cover algorithm to design peptides. The input is the directory where the cluster files are located. Here we will run SW\_SC.py on all files located in the clusters folder. The -e flag can be used to exclude any peptides that contain any non-AA characters such as “X”. To ensure complete coverage of all epitopes of a certain size in the sliding window design, the ideal step size will be: window size - (number of AA in target epitope - 1). Here, our target epitope size is 9. [SW_SC.py]

Command:
```
SW_SC.py -u /clusters/id_70_SW_SC_w30_s22_x9_sumStats.tsv -e x -s 22 -x 9 -y 30  clusters/*
```

Concatenate SW_SC.py output fasta files for each cluster into a single fasta file.

Command:
```
cat *.fasta > poxviridae_id70_all_SWSC-x9y30.fasta
```
Convert fasta file into .csv file with first column as Probe\_id (ex. POX\_000001, POX\_000002…) and the second column as the peptide sequences. The links between coded peptide names and the names in the fasta file will be saved to a map file “POX_map.csv”. This map file can later be used to build out a metadata file for the library.

Command:
```
fasta_to_encodable.py -i poxviridae_id70_all_SWSC-x9y30.fasta -o POX_encodable.csv -p POX -m POX_map.csv
```
Generate nucleotide encodings for all designed peptides [Step 1: [oligo\_encoding], Step 2: oligo_encoding.py]. 

Step 1:

In step 1, a user defined number of possible encodings are randomly generated for each peptide and then the top n encodings are selected depending on the user-specified target GC content.

Here we create 10,000 encodings for each peptide, but only output the 300 encodings that have the lowest absolute deviation from the specified GC content (0.55). We also use two cores for this analysis, as denoted by the '-c 2' option. The input\_file must contain lines of the form {seq},{name}, where the length of each line can be a maximum of 128. After completion, ‘output\_ratio.csv’ will contain the necessary information for input to the neural network in step 2 and ‘out_seqs.csv’ will contain the top 300 encodings for each sequence.
Command:

```
oligo_encoding_ \
-r output_ratio.csv \
-s out_seqs.csv \
-n 300 \
-c 2 \
-p Library-Design/scripts/oligo_encoding/codon_weights_test.csv \
-i POX_encodable.csv \
-t 10000 \
-g 0.55
```
Step 2:
In step 2, the top encodings from step 1 are input into a deep learning model to score and select the top n encodings based on 88 features and trained using relative peptide abundances observed within PepSeq libraries (64 codons, 4 NT, 20AA).

Using the previously created 'output\_ratio' and 'out\_seqs' files, use deeplearning\_model to predict the best encodings in out\_seqs using the data in output\_ratio. In this example, 10 sequences (of k encodings each) will be processed by the Neural Network at a time. The --read\_per_loop flag should be lowered for machines with less memory, and can be increased on machines with more. The encodings with the top 3 lowest absolute neural network predictions will be output for each input sequence.


Command:
```
encoding_with_nn.py \
-m DeepLearning_model_R_1539970074840_1_20181019 \
-r output_ratio \
-s out_seqs \
-o POX_best_encodings.csv \
--subsample 300 \
--read_per_loop 10 \
-n 3
```

The final output will be a comma separated values (csv) with the top 3 encodings for each peptide.

