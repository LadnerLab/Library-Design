# Swiss KB to Fasta -- Translate data files from Uniprot into Fasta Files
## Introduction
This script is used to translate data files that are downloaded from Uniprot into fasta files. 
Different options enable control over the output sequences, including the different tags that will be placed in the name of 
each output sequence. 

## Usage:

```
usage: swisskb_to_fasta.py [-h] [-s SWISS] [-o OUTPUT] [-t TAGS]
                           [-l RANKED_LINEAGE] [-m RANK_MAP]

Script to convert UniprotKB/Swiss-Prot entries to a fasta file. The names for
each entry can be specified from the tags found in each entry.

optional arguments:
  -h, --help            show this help message and exit
  -s SWISS, --swiss SWISS
                        Name of file containing UniprotKB/Swiss-prot entries.
  -o OUTPUT, --output OUTPUT
                        Name of file to write outputs to, data in file will be
                        formatted in the FASTA format.
  -t TAGS, --tags TAGS  Tags that will be included in each sequence name, note
                        that id will always be collected and used for sequence
                        names.To include multiple tags, provide this argument
                        multiple times, each time providing a tag to include.
                        Note that the OXX tag is not part of the uniprot
                        standard, and both OXX/OX may be included. When OXX is
                        included, each entry in the output will include the
                        following information:
                        OXX=original_id,species_id,genus_id,family_id. Note
                        that if OXX is included, both the ranked_lineage and
                        rank_map flags will need to be included. Possible
                        values for tags include: AC, DE, DR, DT, FT, GN, ID
                        (always included), KW, OC, OH, OS, OX, PE, RA, RC, RG,
                        RL, RN, RP, RT, RX, SQ, OXX.
  -l RANKED_LINEAGE, --ranked_lineage RANKED_LINEAGE
                        Map containing taxid|tax_info pairings which can be
                        used to create OC tags. Note that is the OC tag is
                        provided from the command line, then this tag must be
                        provided as well. Inclusion of this argument without
                        the OC tag will be ignored.
  -m RANK_MAP, --rank_map RANK_MAP
                        Map containing taxid|tax_rank pairings, which will be
                        be used to annotate taxid taxonomic rank info. Note
                        that this argument is optional when combined with the
                        'OC' tag. Note that this argument will be ignored if
                        'OC' and '--ranked_lineage' are not also provided. If
                        provided, each OC tag will be formatted
                        'OC=12345,FAMILY', where 'FAMILY' is the taxonomic
                        rank identifier for the id.
```


## Example Command: 

``` console
./swisskb_to_fasta.py 
-t ID -t AC -t OXX 
-l rankedlineage.dmp 
-s uniprot_trembl_viruses.dat 
-o uniprot_trembl_viruses.fasta 
-m taxid_taxrank_map.dmp
```
