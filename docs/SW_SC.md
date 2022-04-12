---
layout: default
title: SW-SC
permalink: /SW_SC/
---
# SW-SC

## Overview
The *SW-SC* script uses a combined sliding window/set cover algorithm to design peptides for a given protein sequence fasta file, or a directory of protein sequence fasta files.

## Inputs

The only required input is a fasta-formatted file containing a set of target protein sequences.

## Outputs

There is one required output, a fasta-formatted file containing the peptide sequences. 

There is one optional output, a tab-delimited summary file, which shows the number of peptides designed for each input file (one line per input file).

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 

### Tutorial/Use

 The input is the directory where the cluster files are located. Here we will run SW_SC.py on all files located in the clusters folder. The -e flag can be used to exclude any peptides that contain any non-AA characters such as “X”. To ensure complete coverage of all epitopes of a certain size in the sliding window design, the ideal step size will be: window size - (number of AA in target epitope - 1). Here, our target epitope size is 9. [SW_SC.py]

Command:
```
SW_SC.py \
  -u  /clusters/id_70_SW_SC_w30_s22_x9_sumStats.tsv \
  -e x -s 22 -x 9 -y 30  clusters/*
```
