#!/usr/bin/env python3

import argparse, os
import kmertools as kt		#Available at https://github.com/jtladner/Modules
import fastatools as ft		#Available at https://github.com/jtladner/Modules
import numpy as np
from collections import defaultdict

#Example command: /Users/colleenung/Documents/GitHub/Library-Design/extensions/design_coverage.py -d /Users/colleenung/OneDrive\ -\ Northern\ Arizona\ University/PepSeq_Designs/PM1/Design/PM1combo_v1v2_rr10-30.fasta -t /Users/colleenung/OneDrive\ -\ Northern\ Arizona\ University/PepSeq_Designs/PM1/Design/targetVerification --swCtoS -k 9

parser = argparse.ArgumentParser(description='A script that will calculate the coverage of given target(s) within the design. Coverage will be calculated on a per sequence basis and overall.')

parser.add_argument("-k", "--ksize", default=9, type=int, metavar='\b', help="Kmer size to use when calculating coverage [default: 9].")
parser.add_argument("-p", "--perseqoutput", default="coverage-per-seq.tsv", metavar='\b', help="Name of output tsv file with calculated kmer coverage for each sequence (> in length than ksize) in input file. [default: coverage-per-seq.tsv]")
parser.add_argument("-s", "--statsoutput", default="coverage-stats.tsv", metavar='\b', help="Name of output tsv file with descriptive statistics for per seq coverage and calculated overall coverage. [default: coverage-stats.tsv]")
parser.add_argument("--swCtoS", default=False, action="store_true", help="Use this flag if Cysteine residues were converted to Serine residues in the SW portion of the design.")
parser.add_argument("-e", "--extensions", default=".fasta,.fna,.ffn,.faa,.frn,.fa", metavar='\b', help="Only target files with these following file extensions will be grabbed. [default: .fasta, .fna, .ffn, .faa, .frn, .fa]")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-d", "--design", required=True, metavar='\b', help="Design file")
reqArgs.add_argument("-t", "--targets", required=True, metavar='\b', help="Target sequences for which coverage will be calculated; Can provide directory, list of paths or singular file path")

args = parser.parse_args()

targetPaths=[]

args.extensions = args.extensions.split(",")

#Determine if targets were provided as a single dir or as a list of paths
if os.path.isdir(args.targets):
	for name in os.listdir(args.targets):
		a= args.targets + "/" + name
		if os.path.isfile(a) and a.endswith(tuple(args.extensions)):
			targetPaths.append(a)
elif "," in args.targets:
	b= args.targets.split(",")
	for a in b:
		targetPaths.append(a.strip())
else:
	targetPaths.append(args.targets)
	
#Reading in metafile and storing each file pair in a dictionary; keys are the path to the original target file, values are the path to file containing designed peptides
# metaDict = defaultdict(str)
# with open(args.metafile, "r") as fin:
# 	for line in fin:
# 		line= line.strip().split("\t")
# 		metaDict[line[0]]= line[1]

#Creating set of all unique kmers within design file
designkSet= kt.kmerSetFasta(args.design, args.ksize, filter=[])
ct=0

#Opening output files
perseqoutputName= str(args.ksize) + "mer-" + args.perseqoutput
statsoutputName= str(args.ksize) + "mer-" + args.statsoutput
fout1= open(perseqoutputName, "w")
fout2= open(statsoutputName, "w")

for targetF in targetPaths:
	targetkmers=[]
	
	#Reading in target fasta file. Returns dictionary containing name:sequence. Values in dict are formatted as a list in case of duplicate names.
	targetseqD= defaultdict(list)
	names, seqs = ft.read_fasta_lists(targetF)
	seqs = [x.upper() for x in seqs]
	c=0
	for n in names:
		targetseqD[n].append(seqs[c])
		c+=1
	
	#Dictionary will store # of unique kmers in each seq
	uniquekmerctD= defaultdict(list)
	
	#Calculating coverage at a sequence level
	coverageperseqD= defaultdict(list)
	for name, s in targetseqD.items():
		for sequence in s:
			if args.swCtoS:
				if "C" in sequence:
					sequence= sequence.replace("C", "S")
			#Creating set of all unique kmers within sequence
			sSet= kt.kmerSet(sequence, args.ksize, filter=["X"])
			
			if len(sSet)>0:
				xmersCovered= sSet.intersection(designkSet)
				percentCovered= (len(xmersCovered) / len(sSet))*100
				coverageperseqD[name].append(float(percentCovered))
				uniquekmerctD[name].append(len(sSet))
			else:
				coverageperseqD[name].append("NA")
				uniquekmerctD[name].append("NA")

	#Calculating overall coverage
	t= kt.kmerSetFasta(targetF, args.ksize, filter=["X"])		#Creating set of all unique kmers within target file
	if args.swCtoS:
		for seq in t:
			if "C" in seq:
				seq= seq.replace("C", "S")
			targetkmers.append(seq)
	targetkSet= set(targetkmers)
	intersect= len(targetkSet.intersection(designkSet))
	overallcoverage= (intersect/len(targetkSet))*100

	#Writing out tsv file with per seq kmer coverage
	if ct == 0:
		header= "Target\tSequence Name\t# of unique %dmers\t%dmer coverage in design\n" % (args.ksize,args.ksize)
		fout1.write(header)
	for name, coverage in coverageperseqD.items():
		d=0
		for c in coverage:
			if c == "NA":
				d+=1
			else:
				uniquekmerct= int(uniquekmerctD[name][d])
				fout1.write("%s\t%s\t%d\t%.3f\n" % (targetF, name, uniquekmerct, c))
				d+=1
	
	#statsoutputName= os.path.splitext(os.path.basename(targetF))[0] + "_coverage-stats.tsv"
	#Preparing coverageperseqD to perform numpy functions
	coverageperseqL=[]
	for coverage in coverageperseqD.values():
		for c in coverage:
			if c != "NA":
				coverageperseqL.append(c)
	
	#Writing out file with descriptive statistics
	if ct == 0:
		header= "Target\tMaximum*\tQ3*\tMedian*\tQ1*\tMinimum*\tIQR*\tMean*\tOverall coverage\t\t*Per seq %dmer coverage" % args.ksize
		fout2.write(header)
	maximum= max(coverageperseqL)
	q3= np.percentile(coverageperseqL, 75, interpolation = 'midpoint')
	median= np.percentile(coverageperseqL, 50, interpolation = 'midpoint')
	q1= np.percentile(coverageperseqL, 25, interpolation = 'midpoint')
	minimum= min(coverageperseqL)
	IQR= q3-q1
	mean= np.mean(coverageperseqL)
	
	line2= "\n%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.3f" % (targetF,maximum,q3,median,q1,minimum,IQR,mean,overallcoverage)		
	fout2.write(line2)

	ct+=1

fout1.close()
fout2.close()
