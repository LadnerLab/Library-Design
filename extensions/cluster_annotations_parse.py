#!/usr/bin/env python3

import argparse, os, re
import fastatools as ft		#Available at https://github.com/jtladner/Modules

#Example command: cluster_annotations_parse.py -f clusters/* -o 10240_Poxviridae_clusterinfo.tsv


parser = argparse.ArgumentParser(description='''A script that will extract annotations from sequence names, given one or more fasta files. 
												Expecting sequences names in the same format as produced by GenBank or UniProt.
Output TSV file will contain the following: cluster #, an empty column for labeling clusters and one column per name annotation. Pair of brackets
following annotation represents the number of sequences with the same annotation in that cluster.''')

parser.add_argument("-o", "--output", default="clusterinfo.tsv", metavar='\b', help="Name of output TSV file containing annotations parsed from input fasta file(s). [default: clusterinfo.tsv]")
#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-f", "--fasta", required=True, nargs="+", metavar='\b', help="Fasta file(s) to parse annotations from. Can take one fasta file or a directory of fasta files.")

args = parser.parse_args()


with open(args.output, "w") as fout:
	fout.write("Cluster#\tProtein\tAnnotations [SeqCount]\n")
	
	for file in args.fasta:
		#Parsing cluster number from file name(s) as UCLUST automatically generates numbered clusters. Assuming cluster number is preceded by underscore (e.g. id_70_24 represents cluster #24).
		searchstr= ".*_(\d+)"
		regexresult= re.search(searchstr, os.path.basename(file))
		clusterNum= regexresult.group(1)
		
		#Reads in fasta file. Returns two lists, the first containing seq names and the second containing its sequences.
		names, seqs = ft.read_fasta_lists(file)
		annotL=[]
		for n in names:
			#Parsing annotations from UniProt seq names
			if n.startswith("tr"):
				annot= re.search("^[^ ]+ (.*) OS=.*", n)
			
			#One format of seqs names from UniProt and GenBank both start with "sp"; differentiating the two and parsing accordingly
			elif n.startswith("sp"):
				#UniProt
				if "OS=" in n:
					annot= re.search("^[^ ]+ (.*) OS=.*", n)
				#GenBank
				elif "RecName:" in n:
					if ";" in n:
						annot= re.search("^.+? Full=(.*?);", n)
					else:
						annot= re.search("^.+? Full=(.*)", n)
				else:
					print("Unable to extract info from %s in cluster %s" % (n, clusterNum))
			
			#Parsing annotations from GenBank seq names, for a variety of formats
			elif n.startswith("pdb"):
				annot= re.search("^[^,]+, (.*)", n)
			elif n.startswith("pir"):
				annot= re.search("^[^ ]+ (.*) -.*", n)
			elif n.startswith("prf"):
				annot= re.search("^[^ ]+ (.*)", n)
			elif n.endswith("]"):
				annot= re.search("^[^ ]+ (.*) \[.*\]", n)		
			
			else:
				print("Unable to extract info from %s in cluster %s" % (n, clusterNum))
				
			if annot is None:
				print("Unable to extract info from %s in cluster %s" % (n, clusterNum))
			if annot is not None:
				annotL.append(annot.group(1))

		seqCount={}
		for a in annotL:
			#If protein name exists in dict, increment value for every occurence. If it doesn't, will create key and assign value of 0 (and then add 1).
			seqCount[a] = seqCount.get(a,0) + 1	
		#Copying the same dictionary, but ordering the annotations by SeqCount (# of sequences with the same annotation) to make manual assignments downstream easier
		seqCountordered={k: v for k, v in sorted(seqCount.items(), reverse=True, key=lambda item: item[1])}
		
		fout.write("%s\t\t" % clusterNum)
		for annotation, counts in seqCountordered.items():
			fout.write("%s [%d]\t" % (annotation, counts))
		fout.write("\n")
