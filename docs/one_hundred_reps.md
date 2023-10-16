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

There is one optional output, a tab-delimited map file, which relates all removed sequences to the representatives to which they are identical (one line per unique representative).

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 
- You will, however, need to visit [this link](https://github.com/LadnerLab/ProteinOligoLibrary/blob/master/protein_oligo_library.py) to download the Protein Oligo Library file and move it to the directory where the scripts you intend to run are located, or place the file in your PYTHONPATH. For example, if you have Anaconda installed, you can move the library to your base environment by putting it here: "~/conda_directory_name/lib/python3.x/site-packages/".

- To install the C version of One Hundred Reps:
    - Download the [source files from GitHub](https://github.com/LadnerLab/Library-Design/tree/master/one_hundred_reps/c)
    - Through the terminal, change your working directory to the location of the source files on your computer.
    - Enter the following command: `make`
    - If successful, this should generate an executable named `one\_hundred\_reps`

## Use

In this example, our [input file](https://github.com/LadnerLab/Library-Design/blob/master/examples/poxviridae_unaligned_min30AA.fasta) will contain a set of 17,786 protein sequences dowloaded from UniProt. 

The [output](https://github.com/LadnerLab/Library-Design/tree/master/examples/expectedOutputs/onehundredreps) should contain 14,505 unique representatives from this input file.

Command (Python version):
```
one_hundred_reps.py \
  -f poxviridae_unaligned_min30AA.fasta \
  -r poxviridae_unaligned_min30AA_100rep.fasta \
  -m 100_rep_map.txt
```

    - This real world example took <2 min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4). 

Command (C version, using 2 threads):
```
one_hundred_reps \
  -f poxviridae_unaligned_min30AA.fasta \
  -n 2 \
  -m 100_rep_map.txt
```

    - This real world example took ~# min to complete on a Macbook Pro laptop (Apple M1, macOS v11.6.4). 
