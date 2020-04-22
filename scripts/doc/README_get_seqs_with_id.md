# Get seqs with ID -- Retrieve the sequences from a FASTA file that contain specified Taxonomic Identifiers 

## Introduction
This script is used to retrieve the sequences from a FASTA file that contain one or more target taxonomic identifiers.
The script uses the ```OXX``` tag that is introduced by [swiss_kb_to_fasta](README_swiss_kb_to_fasta.md).
**Note**: This script does not consider the position of the identifiers, it parses out all of the IDs, 
and any sequence whose ids overlap with those in the target are included in the output. 

## Script location:
[get_seqs_with_id.py](../get_seqs_with_id.py "The script can be found here")


## Usage: 
```
usage: get_seqs_with_id.py [-h] [-i IDS] [-f FASTA]

Retrieve all sequences in a file that have contain a specified taxonomic id.

optional arguments:
  -h, --help            show this help message and exit
  -i IDS, --ids IDS     Name of file containing one id per line
  -f FASTA, --fasta FASTA
                        Name of fasta file whose sequences will be searched
                        for the specified ids.
```

## Example Command: 
``` console
    ./get_seqs_with_id.py -i target_ids.txt
    -f output_fasta.fasta
```

The file ```target_ids.txt``` looks something like:
```
1159901
6399396
```
