#!/usr/bin/env python3

import argparse
import fastatools as ft
import kmertools as kt
from collections import defaultdict
import os
import subprocess

parser = argparse.ArgumentParser(description='''A script that will compare two input fasta files using kmers and generate two combo versions of the fasta files. Outputs will also include two other fasta files, one with sequences from file 1 that contain kmers not present in file 2 and vice versa.''')

parser.add_argument("-k", "--ksize", default=7, type=int, metavar='\b', help="kmer size to use for comparison [default: 7].")
#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-u", "--uniprot", required=True, metavar='\b', help="Fasta file from UniProt")
reqArgs.add_argument("-g", "--genbank", required=True, metavar='\b', help="Fasta file from GenBank")

args = parser.parse_args()


#Reads in fasta file, returns dictionaries containing name:seq and name:kmers
def fastaDict_kmerDict(file):
	fastaDict= ft.read_fasta_dict_upper(file)
	kmerDict= defaultdict(list)

	for name, seq in fastaDict.items():
		for i in range(len(seq)-args.ksize+1):
			kmerDict[name].append(seq[i:i+args.ksize])
	return fastaDict, kmerDict

fastaDict_UP, kmerDict_UP= fastaDict_kmerDict(args.uniprot)
fastaDict_GB, kmerDict_GB= fastaDict_kmerDict(args.genbank)

#Creating an inverse dictionary containing kmers:name
for key, value in kmerDict_UP.items():
	kmerDict_UP[key]= ", ".join(value)
for key, value in kmerDict_GB.items():
	kmerDict_GB[key]= ", ".join(value)

invkmerDict_UP = defaultdict(list)
invkmerDict_GB = defaultdict(list)
{invkmerDict_UP[v].append(k) for k, v in kmerDict_UP.items()}
{invkmerDict_GB[v].append(k) for k, v in kmerDict_GB.items()}
# invkmerDict_UP= {v: k for k, v in kmerDict_UP.items()}
# invkmerDict_GB= {v: k for k, v in kmerDict_GB.items()}


#Set of unique kmers for each file, filtering out kmers that contain "X" which are ambiguous amino acids
uniqkmerSet_UP= kt.kmerSetFasta(args.uniprot, args.ksize, filter=["X"])
uniqkmerSet_GB= kt.kmerSetFasta(args.genbank, args.ksize, filter=["X"])

#Comparing kmers from each file; any sequences with kmers exclusive to one file are added to a dictionary
exclusivetoUP= uniqkmerSet_UP.difference(uniqkmerSet_GB)
exclusivetoGB= uniqkmerSet_GB.difference(uniqkmerSet_UP)

outDict_UP= defaultdict(str)
for x in exclusivetoUP:
	for kmers, name in invkmerDict_UP.items():
		if x in kmers:
			for n in name:
				outDict_UP[n]= fastaDict_UP[n]
outDict_GB= defaultdict(str)
for x in exclusivetoGB:
	for kmers, name in invkmerDict_GB.items():
		if x in kmers:
			for n in name:
				outDict_GB[n]= fastaDict_GB[n]


#Checking if "UP_exclusive" and "GB_exclusive" directories exist; if not, create them
d1= "./UP_exclusive"
d2= "./GB_exclusive"
if not os.path.exists(d1) and not os.path.exists(d2):
	os.mkdir(d1)
	os.mkdir(d2)
	
#Writing out two fasta files to above directories, one with sequences from file 1 that contain kmers not present in file 2 and vice versa
fp1= d1 + "/uniq_" + os.path.basename(args.uniprot)	#Given a path, os.path.basename() will extract the file name and return it as a string
fp2= d2 + "/uniq_" + os.path.basename(args.genbank)	#Given a path, os.path.basename() will extract the file name and return it as a string
ft.write_fasta_dict(outDict_UP, fp1) 
ft.write_fasta_dict(outDict_GB, fp2)


#Checking if the combo directories exist; if not, create them
cd1= "./cmbo_uniqGB_to_UP"
cd2= "./cmbo_uniqUP_to_GB"
if not os.path.exists(cd1) and not os.path.exists(cd2):
	os.mkdir(cd1)
	os.mkdir(cd2)

#Function to append to a fasta file
def append_fasta_dict(fD, filename):
    fout=open(filename, 'a')
    for k,v in fD.items():
        fout.write(">%s\n%s\n" % (k, v))
    fout.close()

#Copying original files to above directories, to which the unique sequences from the other database will be appended to to create the combo files
cfp1= cd1 + "/cmbo_" + os.path.basename(args.uniprot)	#Given a path, os.path.basename() will extract the file name and return it as a string
cfp2= cd2 + "/cmbo_" + os.path.basename(args.genbank)	#Given a path, os.path.basename() will extract the file name and return it as a string
args1= "cp %s %s" % (args.uniprot, cfp1)
args2= "cp %s %s" % (args.genbank, cfp2)
subprocess.call(args1, shell=True)
subprocess.call(args2, shell=True)
append_fasta_dict(outDict_GB, cfp1)
append_fasta_dict(outDict_UP, cfp2)

