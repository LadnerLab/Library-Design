#!/usr/bin/env python3

import fastatools as ft  # Available here: https://github.com/jtladner/Modules
import argparse, random


def main():
    arg_parser = argparse.ArgumentParser( description = "Mutate input sequences to generate diverse datasets.", formatter_class=argparse.ArgumentDefaultsHelpFormatter )

    arg_parser.add_argument( '-n', '--num', help = "Number of mutated sequences to output per input sequence.", default = 30, type = int )

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-i', '--input', help = "Fasta file contianing the protein sequence(s) to mutate.", required=True )
    reqArgs.add_argument( '-o', '--output', help = "Base name for output fasta files.", required=True )
    reqArgs.add_argument( '-d', '--diverg', help = "Level of divergence from the input sequence. Should be between 0 and 1. Can include multiple comma-delimited values.", required=True )

    args = arg_parser.parse_args()
    
    # Possible amino acids
    AAs = ['A', 'C', "D", "E", 'F', "G", "H", 'I', "K", 'L', 'M', 'N', "P", "Q", "R", "S", "T", 'V', 'W', 'Y']
    
    # Parse target divergences
    divergs = [float(d) for d in args.diverg.split(",")]
    
    # Read in input seqs
    iD = ft.read_fasta_dict_upper(args.input)
    
    for d in divergs:
        outN = []
        outS = []
    
        for n,s in iD.items():
            newS = s
            c=0
            muts = int(d*len(s))
            while c<args.num:
                c+=1
                
                sites = random.choices(range(len(s)), k=muts)
                for site in sites:
                    subAAs = AAs[::]
                    subAAs.remove(newS[site])
                    newS = newS[:site] + random.choice(subAAs) + newS[site+1:]
                
                outS.append(newS)
                outN.append("%s_d%.3f_%03d" % (n, d, c))
            
        ft.write_fasta(outN, outS, "%s_d%.3f.fasta" % (args.output, d))
    
if __name__ == '__main__':
    main()
