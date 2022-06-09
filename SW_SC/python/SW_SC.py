#!/usr/bin/env python

import argparse, random, os
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files. Facilitates batch processing. Output names will be generated using each input file name.", nargs="*")
    parser.add_argument( "-u", "--summary", help="Name for a tab-delimited output file summarizing the number of peptides designed for each input set of targets.")
    parser.add_argument("-i", "--inp", help="Input file name. Should contain target protein sequences from which to design peptides. Can be used along with -o if designing for a single target set.")
    parser.add_argument("-o", "--out", help="Output file name. Will be a list of peptides, 1 per line. Can be used along with -i if designing for a single target set.")
    parser.add_argument( '-s', '--step_size', help = "Number of amino acids to move between each window.", default = 1, type = int )
    parser.add_argument("-t", "--target", default=1, type=float, help="Target xmer coverage. Algorithm will continue until at least this proportion of total Xmers are in the design. If '--pre' option is used, Xmers in these predesigned peptides will also be considered in this threshold.")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-x", "--xMerSize", type=int, help="Size of Xmers, which represent potential linear epitopes contained within peptides/Yemrs.", required=True)
    reqArgs.add_argument("-y", "--yMerSize", type=int, help="Size of Ymers, which represent potential peptides for inclusion in the assay.", required=True)

    args = parser.parse_args()
    
    #Create set of characters to exclude
    args.exSet = set(args.exclude)
    
    # Open output summary file for writing, if requested
    if args.summary:
        fout = open(args.summary, "w")
        fout.write("File\tNumPeps\n")
    
    #Run set cover analyses
    for each in args.inputs:
        numPep = design(each, "%s_SWSC-x%d-y%d.fasta" % (os.path.basename(each), args.xMerSize, args.yMerSize), args)

        if args.summary:
            fout.write("%s\t%d\n" % (each, numPep))
    
    if args.inp and args.out:
        numPep = design(args.inp, args.out, args)

        if args.summary:
            fout.write("%s\t%d\n" % (args.inp, numPep))

#----------------------End of main()

def design(inp, out, args):

    # Generate dict with xmer counts
    xcD = {}
    
    tN, tS = ft.read_fasta_lists(inp)
    for s in tS:
        xL = kt.kmerList(s, args.xMerSize)
        for x in xL:
            if len(set(x).intersection(args.exSet)) == 0:
                xcD[x] = xcD.get(x, 0) + 1

    #Save count of total xmers in targets
    totalX = len(xcD)

    # Score each target sequence by summing contained xmer scores
    maxScore = 0
    repS = ""
    repN = ""
    for i,s in enumerate(tS):
        theseXs = kt.kmerList(s, args.xMerSize)
        thisScore = sum([xcD[x] for x in theseXs if x in xcD])
        if thisScore > maxScore:
            maxScore = thisScore
            repS = s
            repN = tN[i]

    # Generate peptides using a sliding window across the chosen representative sequence
    rep = [Sequence( name = repN, sequence = repS )]
    designer = LibraryDesigner( window_size = args.yMerSize, step_size = args.step_size )
    library = designer.design( rep )

    repD = {e.name:e.sequence for e in library}
    repNames = sorted(list(repD.keys()))
    repSeqs = [repD[n] for n in repNames]

    
    # Remove xmers covered by the sliding window peptides
    for s in repSeqs:
        xL = kt.kmerList(s, args.xMerSize)
        for x in xL:
            if x in xcD:
                del(xcD[x])
        
    # Read in all yMers in targets
    ysD = {}
    yNameD = {}
    for i,s in enumerate(tS):
        yL = kt.kmerList(s, args.yMerSize)
        for j, y in enumerate(yL):
            if len(set(y).intersection(args.exSet)) == 0:
                ysD[y] = 0
                yNameD[y] = "%s_%04d" % (tN[i], j)
    
    # Design peptides
    newSeqs = []
    newNames = []
    
    while (1-(len(xcD)/totalX)) < args.target:
        
        thisPep = choosePep(ysD, xcD, args)
        thisName = yNameD[thisPep]
        newSeqs.append(thisPep)
        newNames.append(thisName)
        
        #Remove selected peptide from ysD
        del(ysD[thisPep])
        
        #Remove covered xMers from xcD
        for eachX in kt.kmerList(thisPep, args.xMerSize):
            if eachX in xcD:
                del(xcD[eachX])

    # Write out peptides
    ft.write_fasta(repNames+newNames, repSeqs+newSeqs, out)
    
    return len(repSeqs+newSeqs)

class LibraryDesigner():
    def __init__( self, window_size = 0, step_size = 0 ):
        self.window_size = window_size
        self.step_size   = step_size

    def _get_oligos( self, sequence ):
        xmers = set()

        start = 0
        end = self.window_size

        seq = sequence.sequence

        if len( seq ) < self.window_size and 'X' not in seq:
            xmers.add( seq )
        while end < len( seq ):
            xmer = seq[ start:end ]

            new_name = sequence.name + "_%03d_%03d" % ( start, end )
            if not 'X' in xmer and '-' not in xmer:
                xmers.add( Sequence( name = new_name, sequence = xmer))
                
                # If the remaining sequence is less than the step size
                # Add a peptide covering the end of the seq
                if  0 < len( seq[end:] ) < (self.step_size):
                    xmer = seq[-self.window_size:]
                    new_name = sequence.name + "_%03d_%03d" % ( len(seq)-self.window_size, len(seq) )
                    xmers.add( Sequence( name = new_name, sequence = xmer))
                
            start += self.step_size
            end   = start + self.window_size
        return xmers

    def design( self, sequences ):
        all_oligos = set()

        for seq in sequences:
            oligos = self._get_oligos( seq )
            all_oligos |= oligos
        return all_oligos

class Sequence:
    def __init__( self, name = "", sequence = "" ):
        self.name     = name
        self.sequence = sequence
    def __hash__( self ):
        return hash( self.sequence )
    def __eq__( self, other ):
        return self.sequence == other.sequence
    def __ne__( self, other ):
        return not self.__eq__( other )
    def __str__( self ):
        return '>%s\n%s\n' % ( self.name, self.sequence )
    def __len__( self ):
        return len( self.sequence )


def choosePep(ysD, xcD, args):
    
#    print(len(ysD), len(xcD))
#    if len(xcD) == 15:
#        print (xcD)
    
    #Calculate scores for xMers
    for y in ysD:
        theseXs = kt.kmerList(y, args.xMerSize)
        ysD[y] = sum([xcD[x] for x in theseXs if x in xcD])

    #Dict by score
    scoreD = defaultdict(list)
    for k,v in ysD.items():
        scoreD[v].append(k)
    
    #Choose peptide
    thisMax = max(scoreD.keys())
    
    thisChoice = random.choice(scoreD[thisMax])
    
    return thisChoice

def writeXmerDict(xD, outname):
    with open(outname, "w") as fout:
        fout.write("Xmer\tCount\n")
        for k,v in xD.items():
            fout.write("%s\t%d\n" % (k, v))

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

