#!/usr/bin/env python3

# By Jason Ladner

# Wrapper for translating a list of genbank assembly accessions into host metadata

import argparse
import inout as io
#from subprocess import Popen, PIPE

def main():

    #To parse command line
    p = argparse.ArgumentParser()
    
    p.add_argument('-i', '--input', help='Input file containing a list of genbank accessions')
    p.add_argument('-o', '--out', help='Name for output file')
    p.add_argument('-m', '--max', default=100, type=int, help='Max number of accessions per query.')
    p.add_argument('-d', '--db', default="protein", help='Database from which accessions were obtained. Current supported databases: protein, assembly.')

    args = p.parse_args()

    accL = io.fileList(args.input, header=False)
    
#    cmd="echo $PATH"
#    print (cmd)
#    catch = Popen(cmd, shell=True)

    print("#!/bin/zsh")
    
    for start in range(0, len(accL), 100):
        queries = ",".join(accL[start:start+args.max])

        if args.db.upper() == "PROTEIN":
            cmd='epost -db protein -format acc -id %s | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion,TaxId,Organism,SubType,SubName >>! %s ' % (queries, args.out)
        elif args.db.upper() == "ASSEMBLY":
            cmd='epost -db assembly -format acc -id %s | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -element Accession,Taxonomy,Organism -group Attribute -if Attribute@harmonized_name -equals "host" -element Attribute >>! %s ' % (queries, args.out)
        else:
            print(args.db, "is not a supported database type.")
        
        print (cmd)
#        catch = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, executable="/bin/zsh")

###-----------------End of main()--------------------------->>>

###------------->>>

if __name__ == "__main__":
    main()
