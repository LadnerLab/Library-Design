# Lib Naive -- Using a Sliding Window to Design a Peptide Library

## Introduction
This script is used to design a set of oligos from a fasta file containing 
aligned sequences. 
This script offers 2 modes:
1. Gap-spanning mode: Any oligos that would contain gaps are discarded.
2. Non-gap-spanning mode: When designing oligos, gaps are skipped. 
The sliding window is extended until the designed oligo contains the target number of AA.

## Script location:
[lib_naive.py](../lib_naive.py "The script can be found here")

## Usage: 
```
usage: lib_naive.py [-h] [-q QUERY] [-o OUTPUT] [-g] [-w WINDOW_SIZE]
                    [-s STEP_SIZE]

Simple creation of a library, both gap spanning and non-gap spanning
algorithms

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        Fasta query file.
  -o OUTPUT, --output OUTPUT
                        Fasta file to output
  -g, --gap_span        Fasta query file.
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Window Size to use for grabbing oligos.
  -s STEP_SIZE, --step_size STEP_SIZE
                        Number of amino acids to step after each window.
```


## Example Command:
```
./lib_naive.py -q ORF7b_aligned.fasta -w 30 -s 22 -g -o output_file
```
