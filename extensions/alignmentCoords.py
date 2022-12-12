#!/usr/bin/env python

import argparse, os
import fastatools as ft	#Available at https://github.com/jtladner/Modules
from collections import defaultdict

#Example command: alignmentCoords.py -a hepacivirus_speciesReps_ginsi_sl.faa -c hepacivirus_unaligned_coords.tsv -t 15 -p hepacivirus_

parser = argparse.ArgumentParser(description='''A script that will generate coordinates for an alignment file, given the alignment file and coordinates of the unaligned seq(s). An output file will be produced for each protein of interest (PROT_coords_alignment.tsv).
		Any missing coordinates for a seq will be substituted with earliest and latest known coord for that protein based on other seqs in the alignment. ''')

parser.add_argument("-t", "--threshold", default=15, type=int, metavar='\b', help="Amount of overlap allowed between proteins. Warning is printed to the screen if coordinates exceed threshold. [default: 15]")
parser.add_argument("-p", "--prefix", help="File prefix for output files, if desired")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-a", "--alignment", required=True, metavar='\b', help="Alignment file to generate coordinates for")
reqArgs.add_argument("-c", "--coords", required=True, metavar='\b', help="Input tab-delimited file with the unaligned seqs' **known** coordinates for proteins of interest. Expecting the following header: File	SeqName	Protein	Start	Stop	RevComp")

args = parser.parse_args()


#Reading in alignment file. Returns dictionary containing sequence name:aligned sequence
alignedSeqs= ft.read_fasta_dict_upper(args.alignment)

annotCoords= defaultdict(list) #Dictionary with known protein coordinates for unaligned sequences
outCoords= defaultdict(list) #Cumulative dictionary that will contain coords for all proteins of interest for all seqs { SeqName : [(Protein,start,stop), (Protein,start,stop), (Protein,start,stop)] }

#Reading in coords file with known coordinates for unaligned seqs.
with open(args.coords, "r") as fin:
	next(fin)		#Skips header line
	for line in fin:
		l= line.rstrip("\n").split("\t")
		SeqName=l[1]
		
		if len(l[3].strip())>0 and len(l[4].strip())>0:
			start=int(l[3])
			stop=int(l[4])
			cnt=0
			outCoordsValue=[]
			#Determining start/stop coords for alignment and adding to dictionary annotCoords (proteinName:coordList)
			for i, c in enumerate(alignedSeqs[SeqName], 1):
				if c != "-":
					cnt+=1
					if cnt == start:
						annotCoords[l[2]].append(i)
						outCoordsValue.append(l[2])
						outCoordsValue.append(i)
					elif cnt == stop:
						annotCoords[l[2]].append(i)
						outCoordsValue.append(i)
						outCoords[SeqName].append(tuple(outCoordsValue))
						break
		else:
			print("error: File %s has one or more incomplete pairs of start/stop coordinates." % os.path.basename(args.coords))
			quit()

#For each protein of interest, selecting the earliest and latest alignment coords to use as the start/stop coords for the sequences w/o known coords
earliestlatestDict={}
for proteinName, coordList in annotCoords.items():
	tple= (min(coordList),max(coordList))
	earliestlatestDict[proteinName]= tple
#Identifying which seqs in the alignment file have missing coords and filling in gaps with earliestlatestDict
##First, looking at seqs that are missing coords for all proteins of interest
protsofInterestL= earliestlatestDict.keys()
for SeqName in alignedSeqs.keys():
	if SeqName not in outCoords.keys():
		for protein in protsofInterestL:
			start= earliestlatestDict[protein][0]
			stop= earliestlatestDict[protein][1]
			val= (protein, start, stop)
			outCoords[SeqName].append(val)
##Then, filling in gaps for seqs that have coords for some, but not all, proteins of interest
for SeqName, tupleList in outCoords.items():
	coveredprotsL=[]
	for tple in tupleList:
		coveredprotsL.append(tple[0])
	for protein in protsofInterestL:
		if protein not in coveredprotsL:
			start= earliestlatestDict[protein][0]
			stop= earliestlatestDict[protein][1]
			val= (protein, start, stop)
			outCoords[SeqName].append(val)


#Writing out coords for alignment file into protein-specific coords files.
headerList=[]
for SeqName, tupleList in outCoords.items():
	for tple in tupleList:
		protein=tple[0]
		start=int(tple[1])
		stop=int(tple[2])
		if args.prefix:
			outputfileName= "%s%s_coords.tsv" % (args.prefix, protein)
		else:
			outputfileName= "%s_coords.tsv" % (protein)
			
		with open(outputfileName, "a") as fout:
			if protein not in headerList:
				fout.write("File\tSeqName\tStart\tStop\tRevComp\n")
				headerList.append(protein)
			fout.write("%s\t%s\t%d\t%d\t0\n" % (args.alignment, SeqName, start, stop))


#Print warning to the screen if protein coord overlap > threshold
exceedThresh= defaultdict(list)
for SeqName, tupleList in outCoords.items():
	orderedbyStart= sorted(tupleList, key=lambda tup: tup[1])	#Reordering proteins by start coord
	for x in range(1, len(tupleList)):
		tuple1= tupleList[x-1]
		tuple2= tupleList[x]
		overlap= tuple1[2]-tuple2[1]
		
		#If overlap exceeds threshold, check how many chars from the alignment are not dashes. If # of chars > threshold, print warning
		if overlap > args.threshold:
			overlapStr= tuple2[tuple2[1]:tuple1[2]]
			overlapCoords=(tuple2[1],tuple1[2])
			cnt=0
			for c in overlapStr:
				if c != "-":
					cnt+=1
			if cnt > args.threshold:
				exceedThreshValue=[tuple1[0], tuple2[0], overlapCoords]
				exceedThresh[SeqName].append(tuple(exceedThreshValue))
if len(exceedThresh) > 0:
	print("Amount of overlap between proteins exceeds threshold (>%d) for the following:\n" % args.threshold)
	for SeqName, tupleList in exceedThresh.items():
		for t in tupleList:
			print("%s\t%s\t%s\t%s\n" % (SeqName, t[0], t[1], t[2],))
