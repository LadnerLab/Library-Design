#!/usr/bin/env python

import argparse, random, os, glob
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files. Facilitates batch processing. Output names will be generated using each input file name.", nargs="*")
    parser.add_argument("-u", "--summary", help="Name for a tab-delimited output file summarizing the number of peptides designed for each input set of targets.")
    parser.add_argument("-s", "--step_size", help = "Number of amino acids to move between each window.", default = 1, type = int )
    parser.add_argument("-t", "--targets", default="0.5,0.75,1", help="Target thresholds for xmer coverage (comma-separated). Algorithm will continue until at least the max proportion of total Xmers are in the design, but will write the design out once each threshold is met.")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-x", "--xMerSize", type=int, help="Size of Xmers, which represent potential linear epitopes contained within peptides/Yemrs.", required=True)
    reqArgs.add_argument("-y", "--yMerSize", type=int, help="Size of Ymers, which represent potential peptides for inclusion in the assay.", required=True)

    args = parser.parse_args()
    
    #Parse target thresholds
    targetThresh = sorted(list(set([float(x) for x in args.targets.split(",")])))
    maxThresh = targetThresh[-1]
    otherThresh = targetThresh[:-1]
    
    #Generate directories for each threshold
    for thr in targetThresh:
        if not os.path.isdir("t%.3f" % (thr)):
            os.mkdir("t%.3f" % (thr))
    
    #Create set of characters to exclude
    args.exSet = set(args.exclude)
    
    # Open output summary file for writing, if requested
    if args.summary:
        fout = open(args.summary, "w")
        fout.write("File\tXmerThreshold\tNumPeps\n")
    
    #Run set cover analyses
    for each in args.inputs:   #Step through each input file
        #Run the design
        numPepsD = design(each, maxThresh, otherThresh[::], args)

        if args.summary:
            for k,v in numPepsD.items():
                fout.write("%s\t%.3f\t%d\n" % (each, k, v))
                
    # Generate concatenated output files for each threshold
    for thr in targetThresh:
        ft.combine_fastafiles(glob.glob("t%.3f/*fasta" % (thr)), "SWSC-x%d-y%d-t%.3f.fasta" % (args.xMerSize, args.yMerSize, thr))
        


#----------------------End of main()

def design(inp, maxThresh, otherThresh, args):

    # Dictionary that will be used to keep track of the number of peptides in each design
    numPepD = {}

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

    # Score each target sequence by summing contained xmer scores. This is to choose the representative for the sliding window portion of the design
    maxScore = 0
    repS = ""
    repN = ""
    for i,s in enumerate(tS):   # Stepping through each target sequence
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
    
    while (1-(len(xcD)/totalX)) < maxThresh:
        
        if len(otherThresh) > 0:
            if (1-(len(xcD)/totalX)) >= otherThresh[0]:
                # Write out peptides for "this" thresh design
                ft.write_fasta(repNames+newNames, repSeqs+newSeqs, "t%.3f/%s_SWSC-x%d-y%d-t%.3f.fasta" % (otherThresh[0], os.path.basename(inp), args.xMerSize, args.yMerSize, otherThresh[0]))
                #Add peptide count to dictionary
                numPepD[otherThresh[0]] = len(repSeqs+newSeqs)
                #Delete this threshold from otherThresh list
                del(otherThresh[0])
        
        thisPep = choosePep(ysD, xcD, args)
        
        if thisPep:
            thisName = yNameD[thisPep]
            newSeqs.append(thisPep)
            newNames.append(thisName)
        
            #Remove selected peptide from ysD
            del(ysD[thisPep])
        
            #Remove covered xMers from xcD
            for eachX in kt.kmerList(thisPep, args.xMerSize):
                if eachX in xcD:
                    del(xcD[eachX])
        
        else:
            print("Unable to cover %d Xmers for %s" % (len(xcD), os.path.basename(inp)))
            xcD={}
        
    # Write out peptides for maxThresh design
    ft.write_fasta(repNames+newNames, repSeqs+newSeqs, "t%.3f/%s_SWSC-x%d-y%d-t%.3f.fasta" % (maxThresh, os.path.basename(inp), args.xMerSize, args.yMerSize, maxThresh))
    numPepD[maxThresh] = len(repSeqs+newSeqs)
    
    # Write out peptides for any remaining other thresholds (will happen, for example, if there is just one target seq and therefore, the SW portion covers all Xmers)
    for every in otherThresh:
        ft.write_fasta(repNames+newNames, repSeqs+newSeqs, "t%.3f/%s_SWSC-x%d-y%d-t%.3f.fasta" % (every, os.path.basename(inp), args.xMerSize, args.yMerSize, every))
        numPepD[every] = len(repSeqs+newSeqs)
    
    return numPepD

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
    try:
        thisMax = max(scoreD.keys())
        if thisMax == 0:
            thisChoice=None
        else:
            thisChoice = random.choice(scoreD[thisMax])
    
    except:
        thisChoice=None
    
    return thisChoice

def writeXmerDict(xD, outname):
    with open(outname, "w") as fout:
        fout.write("Xmer\tCount\n")
        for k,v in xD.items():
            fout.write("%s\t%d\n" % (k, v))

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

