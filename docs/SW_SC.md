---
layout: default
title: SW-SC
permalink: /SW_SC/
---
# SW-SC

## Overview
For each cluster, use a combined sliding window/set cover algorithm to design peptides. The input is the directory where the cluster files are located. Here we will run SW_SC.py on all files located in the clusters folder. The -e flag can be used to exclude any peptides that contain any non-AA characters such as “X”. To ensure complete coverage of all epitopes of a certain size in the sliding window design, the ideal step size will be: window size - (number of AA in target epitope - 1). Here, our target epitope size is 9. [SW_SC.py]

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 

### Tutorial/Use
Command:
```
SW_SC.py \
  -u  /clusters/id_70_SW_SC_w30_s22_x9_sumStats.tsv \
  -e x -s 22 -x 9 -y 30  clusters/*
```
Concatenate SW_SC.py output fasta files for each cluster into a single fasta file.

Command:
```
cat *.fasta > poxviridae_id70_all_SWSC-x9y30.fasta
```

Convert fasta file into .csv file with first column as Probe_id (ex. POX_000001, POX_000002…) and the second column as the peptide sequences. The links between coded peptide names and the names in the fasta file will be saved to a map file “POX_map.csv”. This map file can later be used to build out a metadata file for the library.

Command:
```
fasta_to_encodable.py -i poxviridae_id70_all_SWSC-x9y30.fasta \
  -o POX_encodable.csv -p POX -m POX_map.csv
```

## Installation
