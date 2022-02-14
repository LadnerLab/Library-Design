# Kmer-based Oligonucleotide Design
Linux/OSX program to design a library of oligos from input sequences. 
Efficient design of oligonucleotide pools that seeks to minimize the number of oligos produced
while maximizing representation of the input sequences.

**Note** that for larger datasets, it is recommended to use this on a linux-based HPC cluster


### Dependencies
 -  Linux/OSX Operating System, all other files/libraries are included
 - GCC

### Input
    1. An unaligned fasta file containing protein sequences
    2. integer ymer-size, number of amino acids each output oligo will be, note that this must be bigger than xmer-size.
    3. integer xmer-size size of xmers to be represented within the ymers

### Output
	1. Fasta file containing the designed library, each sequence name is the original name of the sequence
       followed by its relative in the input sequence.

### Usage
To get usage info:
```
./set_cover -h
```

#### Example
Create a design of 19-mers, with maximization of the representation of 10-mers within those 19-mers, from seqs.fasta. 
Use 4 threads to speed up the design, and write the best of two designs to oligos.fasta
```
./set_cover -q seqs.fasta -o oligos -x 10 -y 19 -i 2 -t 4
```

### Installation
#### MacOS
You must have make and gcc installed in order to build this program,
this can be achieved with 
```
xcode-select --install
```
from the command line, or through Xcode itself. Once you have these installed,
follow the instructions for Linux.
#### Linux
```
git clone https://github.com/LadnerLab/Library-Design.git
cd setCover/c/
make optimized
```
This produces the set_cover executable, which can be used as described above, with the options described below.

### Options
```

Usage: ./set_cover [ options ]
 -h, --help                 display this help and exit.
 -x                         integer xmer window size. [None, Required]

 -y                         integer ymer window size. [None, Required]

 -e                         A fasta file containing previously
                            designed peptides. Xmers contained in these sequences
                            will not contribute to Ymer scoring in design.

 -r                         default redundancy of xmers in the table. This value grows by one each time a 
 	                        certain xmer is present in the design [1]

 -i                         number of times to do the design, 
 	                        the best one is picked and written to the output files. [1]

 -q                         fasta query file to perform operations on. [None, Required]. 

 -o                         name of file to output to [output.fasta]

 -t                         number of threads to use [1]

 -p                         include this flag in order to perform permutation of xmer functional groups

 -c                         floating point minimum xmer coverage [1]

 -b                         blosum matrix to be used in inclusion of xmer functional groups.
                            Note that blosum90 and blosum62 are hard-coded into this program,
                            and are specified by blosum90 or blosum62. Otherwise, specify a 
                            text file containing a blosum matrix.

 -n                         integer cutoff for whether an amino acid can be substituted 
                            only relationships greater to or equal to this number will be added [0] 

```
