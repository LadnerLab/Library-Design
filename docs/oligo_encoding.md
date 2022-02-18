---
layout: default
title: Oligo Encoding
permalink: /oligo_ecoding/
---
# Oligo Encoding

## Overview
Generate nucleotide encoding for all designed peptides [Step 1: [oligo_ecoding], Step 2: oligo_encoding.py].

### Tutorial/Use

Step 1:
In step 1, a user defined number of possible encodings are randomly generated for each peptide and then the top n encodings are selected depending on the user-specified target GG content.

Command:
```
main \
    -r output_ratio.csv \
    -s out_seq.csv \
    -n 300 \
    -c 2
    -p Library-Design/scripts/oligo_encoding/codon_weights_test.csv \
    -i POX_encodable.csv \
    -t 10000
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
  --subsample 300
  --read_per_loop 10 \
  -n 3
```

## Installation
