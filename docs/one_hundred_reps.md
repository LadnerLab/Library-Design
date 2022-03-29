---
layout: default
title: One Hundred Reps
permalink: /one_hundred_reps/
---

# One Hundred Reps

## Overview

The *One Hundred Reps* scripts can be used to reduce the size of your Target protein sequence file by removing identical protein sequences. This includes perfect matches as well as those that are a subset of another target. A single, longest representative will be maintained in the output.

Two versions of this script are available: a [Python version](https://github.com/LadnerLab/Library-Design/tree/master/one_hundred_reps/python) and a [C version](https://github.com/LadnerLab/Library-Design/tree/master/one_hundred_reps/c). The C version is recommended for large datasets.

## Inputs

The only required input is a fasta-formatted file containing a set of target protein sequences.

## Outputs

There is one required output, a fasta-formatted file containing the unique, representative sequences from the input. 

There is one optional output, a tab-delimited map file, which relates all removed sequences to the representatives to which they are identical. 

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 

- To install the C version of One Hundred Reps:
    - Download the [source files from GitHub](https://github.com/LadnerLab/Library-Design/tree/master/one_hundred_reps/c)
    - Through the terminal, change your working directory to the location of the source files on your computer.
    - Enter the following command: `make`
    - If successful, this should generate an executable names "one\_hundred\_reps"

## Tutorial/Use



Command:
```
onehundredreps.py \
  -f poxviridae_unaligned_30AA.fasta \
  -r poxviridae_unaligned_30AA_100rep.fasta \
  -m 100_rep_map.txt
```

