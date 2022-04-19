---
layout: default
title: Oligo Encoding
permalink: /oligo_encoding/
---
# Oligo Encoding
## Overview

The *oligo_encoding* scripts can be used generate nucleotide encoding for all designed peptides [Step 1: oligo_encoding, Step 2: encoding_with_nn.py].

This script is split into two scripts: Step 1 a [Python version](https://github.com/LadnerLab/Library-Design/blob/master/oligo_encoding/encoding_with_nn.py) and Step 2 a [C version](https://github.com/LadnerLab/Library-Design/tree/master/oligo_encoding). The C version is recommended for large datasets.

## Inputs

The only required input is a fasta-formatted file containing a set of peptides of equal length.

## Outputs

There is one required output, a fasta-formatted file containing the unique, representative sequences from the input. 

There is one optional output, a tab-delimited map file, which relates all removed sequences to the representatives to which they are identical (one line per unique representative).

## Installation

- Because Python is an interpreted language, there is no installation required for the Python version of this script. The only requirement is Python 3. 

- To install the C version of One Hundred Reps:
    - Download the [source files from GitHub](https://github.com/LadnerLab/Library-Design/tree/master/oligo_encoding)
    - Through the terminal, change your working directory to the location of the source files on your computer.
    - Enter the following command: `make`
    - If successful, this should generate an executable named `oligo_encoding`

### Use

Step 1:
In step 1 of this example, a user defined number of possible encodings are randomly generated for each peptide and then the top n encodings are selected depending on the user-specified target GG content.

Command:
```
oligo_encoding \
    -r output_ratio.csv \
    -s out_seq.csv \
    -n 300 \
    -c 2 \
    -p Library-Design/scripts/oligo_encoding/codon_weights_test.csv \
    -i POX_encodable.csv \
    -t 10000 \
    -g 0.55
```

Step 2, the top encodings from step 1 are input into a deep learning model to score and select the top n encodings based on 88 features and trained using relative peptide abundances observed within PepSeq libraries (64 codons, 4 NT, 20AA).

Using the previously created 'output_ratio' and 'out_seqs' files, use deeplearning_model to predict the best encodings in out_seqs using the data in output_ratio. In this example, 10 sequences (of k encodings each) will be processed by the Neural Network at a time. The --read_per_loop flag should be lowered for machines with less memory, and can be increased on machines with more. The encodings with the top 3 lowest absolute neural network predictions will be output for each input sequence.

Command:
```
encoding_with_nn.py \
  -m DeepLearning_model_R_1539970074840_1_20181019
  -r output_ratio \
  -s out_seqs \
  -o POX_best_encoding.csv \
  --subsample 300 \
  --read_per_loop 10 \
  -n 3
```