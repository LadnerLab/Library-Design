---
layout: default
title: Set Cover
permalink: /set_cover/
---
# Set Cover

## Overview

The *Set Cover* scripts can be used to run a set cover algorithm to design peptides for a given protein sequence fasta file, or a directory of protein sequence fasta files.

Two versions of this script are available: a [Python version](https://github.com/LadnerLab/Library-Design/blob/master/setCover/python) and a [C version](https://github.com/LadnerLab/Library-Design/tree/master/setCover/c). The C version is recommended for large datasets.

## Inputs

The only required input is a fasta-formatted file containing a set of target protein sequences.

## Outputs

There is one required output, a fasta-formatted file containing the peptide sequences. 

There is one optional output, a tab-delimited summary file, which shows the number of peptides designed for each input file (one line per input file).

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 

- To install the C version of One Hundred Reps:
    - Download the [source files from GitHub](https://github.com/LadnerLab/Library-Design/tree/master/setCover/c)
    - Through the terminal, change your working directory to the location of the source files on your computer.
    - Enter the following command: `make`
    - If successful, this should generate an executable named `setCover`

## Use

In this example, our [input directory](https://github.com/LadnerLab/Library-Design/tree/master/examples/clusters) will contain ten files of aligned clusters created from a downsampled set of proteins from Poxviridae. (see full tutorial for creation of these files)

The [output](https://github.com/LadnerLab/Library-Design/tree/master/examples/expectedOutputs/setCover) should contain ten files with 30 amino acid long peptides covering the input sequences along with a summary file showing the number of peptides designed for each cluster file. The ten peptide containing fasta files can then be concatenated into a single fasta file.

Command (Python version):
```
setCover.py \
-u poxviridae_id_70_SC_x9_y30_sumStats.tsv \
-x 9 \
-y 30 \
clusters/POX* 
```

    - This real world example took <2 min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4). 

Command (C version, using 2 threads):
```
Coming Soon!
```

    - This real world example took ~# min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4). 
    
The resulting peptide containing fasta files will need to be concatenated into one file containing all of the peptides.

Example Command (Linux):
```cat *.fasta > poxviridae_id70_all_SC-x9y30.fasta 
```