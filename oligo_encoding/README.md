# Oligonucleotide Encoding of Peptides
Scripts to encode a library of peptides as oligonucleotides for PepSeq library synthesis. 

## GPL-3.0-or-later

Visit our [GitHub Pages website](https://ladnerlab.github.io/Library-Design/oligo_ecoding/)

### System requirements
- Linux

#### Tested versions
- Linux: Mint 19.3, CentOS 7

### Software dependencies
    - OpenMP 
    - g++
    - Python 3.5 or greater
    - Pandas python module
    - If you plan to use the provided DeepLearning model, h2o 3.20.0.8 MUST be installed on your system.

### Installation
In order to generate encodings for your protein sequences, the Cpp files must
be compiled. To do this:
```
    make optimized
```
This produces the executable called 'oligo_encoding'.

- Typical installation time <5 minutes. 

### Tutorials

- See [here](https://ladnerlab.github.io/Library-Design/oligo_ecoding/) for SW-SC tutorials.
