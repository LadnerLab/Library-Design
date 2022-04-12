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

## Tutorial/Use

In this tutorial, our [input file](https://github.com/LadnerLab/Library-Design/blob/master/examples/poxviridae_unaligned_min30AA.fasta) will contain a set of 17,786 protein sequences dowloaded from UniProt. 

The [output](https://github.com/LadnerLab/Library-Design/blob/master/examples/expectedOutputs/onehundredreps/poxviridae_unaligned_min30AA_100rep.fasta) should contain 14,505 unique representatives from this input file.

Command:
```SW_SC.py
slidingWindow.py \
-u ./clusters/id_70_SW_w30_s22_sumStats.tsv \
-o poxviridae_id70_all_SW-w30s22.fasta \
-w 30 \
-s 22 \
clusters/* 
```

    - This real world example took <2 min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4).

The resulting peptide containing fasta files will need to be concatenated into one file containing all of the peptides.

Example Command (Linux):
```cat *.fasta > poxviridae_id70_all_SW-w30s22.fasta 
```