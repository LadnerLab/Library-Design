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

## Tutorial/Use

In this tutorial, our [input file](https://github.com/LadnerLab/Library-Design/blob/master/examples/poxviridae_unaligned_min30AA.fasta) will contain a set of 17,786 protein sequences dowloaded from UniProt. 

The [output](https://github.com/LadnerLab/Library-Design/blob/master/examples/expectedOutputs/onehundredreps/poxviridae_unaligned_min30AA_100rep.fasta) should contain 14,505 unique representatives from this input file.

Command (Python version):
```SW_SC.py \
-u ./clusters/id_70_SW_w30_s22_sumStats.tsv \
-o poxviridae_id70_all_SW-w30s22.fasta \
-x 9 \
-y 30 \
clusters/* 
```

    - This real world example took <2 min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4). 

Command (C version, using 2 threads):
```

```

    - This real world example took ~# min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4). 
