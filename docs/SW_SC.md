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

### Use

In this example, our [input directory](https://github.com/LadnerLab/Library-Design/tree/master/examples/clusters) will contain ten files of aligned clusters created from a downsampled set of proteins from Poxviridae. (see full tutorial for creation of these files)

The [output](https://github.com/LadnerLab/Library-Design/tree/master/examples/expectedOutputs/SW_SC)should contain ten files with 30 amino acid long peptides tiling across the input sequences along with a summary file showing the number of peptides designed for each cluster file. The ten peptide containing fasta files can then be concatenated into a single fasta file.

Command:
```
SW_SC.py \
  -u  /clusters/id_70_SW_SC_w30_s22_x9_sumStats.tsv \
  -e x \
  -s 22 \
  -x 9 \
  -y 30  \
  clusters/*
```
