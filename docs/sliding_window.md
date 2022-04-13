---
layout: default
title: Sliding Window
permalink: /sliding_window/
---
# Sliding Window

## Overview

The *Sliding Window* script can be used to tile across input proteins sequences to create peptides of user defined size. The coverage can be set by the user by adjusting the step size between windows.

A python version of this script is available: a [Python version](https://github.com/LadnerLab/Library-Design/blob/master/slidingWindow/python/slidingWindow.py).

## Inputs

The only required input is a fasta-formatted file containing a set of target protein sequences.

## Outputs

There is one required output, a fasta-formatted file containing the peptide sequences. 

There is one optional output, a tab-delimited summary file, which shows the number of peptides designed for each input file (one line per input file).

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 

## Use

In this example, our [input directory](https://github.com/LadnerLab/Library-Design/tree/master/examples/clusters/aligned) will contain ten files of aligned clusters created from a downsampled set of proteins from Poxviridae. (see full tutorial for creation of these files)

The [output](https://github.com/LadnerLab/Library-Design/tree/master/examples/expectedOutputs/slidingWindow)should contain ten files with 30 amino acid long peptides tiling across the input sequences along with a summary file showing the number of peptides designed for each cluster file. The ten peptide containing fasta files can then be concatenated into a single fasta file.

Command:
```
slidingWindow.py \
-u poxviridae_id_70_SW_w30_s22_sumStats.tsv \
-w 30 \
-s 22 \
clusters/POX* 
```

    - This real world example took <2 min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4).

The resulting peptide containing fasta files will need to be concatenated into one file containing all of the peptides.

Example Command (Linux):
```cat *.fasta > poxviridae_id70_all_SW-w30s22.fasta 
```