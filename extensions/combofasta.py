#!/usr/bin/env python3

import argparse, os, subprocess
import fastatools as ft	#Available at https://github.com/jtladner/Modules
import kmertools as kt	#Available at https://github.com/jtladner/Modules
from collections import defaultdict

parser = argparse.ArgumentParser(description='''A script that will compare two input fasta files using kmers and generate two combo versions of the fasta files. Outputs will also include two other fasta files, one with sequences from file 1 that contain kmers not present in file 2 and vice versa.''')

parser.add_argument("-k", "--ksize", default=9, type=int, metavar='\b', help="kmer size to use for comparison [default: 9].")
parser.add_argument("-s", "--subsetMode", default=False, action="store_true", help="This flag can be used if you only wish to check if file2 is a subset of file1. In this mode, only one combo fasta file version will be created— seqs from file1 + seqs w kmers unique to file2. A seperate fasta file containing only the latter will be created as well (file2 seqs that contain kmers not present in file1). [default: None]")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-f1", "--file1", required=True, metavar='\b', help="Fasta file #1")
reqArgs.add_argument("-f2", "--file2", required=True, metavar='\b', help="Fasta file #2")

args = parser.parse_args()


#Reads in fasta file, returns dictionaries containing name:seq and name:kmers
def fastaDict_kmerDict(file):
	fastaDict= ft.read_fasta_dict_upper(file)
	kmerDict= defaultdict(list)

	for name, seq in fastaDict.items():
		for i in range(len(seq)-args.ksize+1):
			kmerDict[name].append(seq[i:i+args.ksize])
	return fastaDict, kmerDict

fastaDict_file1, kmerDict_file1= fastaDict_kmerDict(args.file1)
fastaDict_file2, kmerDict_file2= fastaDict_kmerDict(args.file2)

#If "--subsetMode" flag is NOT provided
if not args.subsetMode:
	#Creating an inverse dictionary containing kmers:name
	for key, value in kmerDict_file1.items():
		kmerDict_file1[key]= ", ".join(value)
	for key, value in kmerDict_file2.items():
		kmerDict_file2[key]= ", ".join(value)

	invkmerDict_file1 = defaultdict(list)
	invkmerDict_file2 = defaultdict(list)
	{invkmerDict_file1[v].append(k) for k, v in kmerDict_file1.items()}
	{invkmerDict_file2[v].append(k) for k, v in kmerDict_file2.items()}
	# invkmerDict_file1= {v: k for k, v in kmerDict_file1.items()}
	# invkmerDict_file2= {v: k for k, v in kmerDict_file2.items()}


	#Set of unique kmers for each file, filtering out kmers that contain "X" which are ambiguous amino acids
	uniqkmerSet_file1= kt.kmerSetFasta(args.file1, args.ksize, filter=["X"])
	uniqkmerSet_file2= kt.kmerSetFasta(args.file2, args.ksize, filter=["X"])

	#Comparing kmers from each file; any sequences with kmers exclusive to one file are added to a dictionary
	exclusivetofile1= uniqkmerSet_file1.difference(uniqkmerSet_file2)
	exclusivetofile2= uniqkmerSet_file2.difference(uniqkmerSet_file1)

	outDict_file1= defaultdict(str)
	for x in exclusivetofile1:
		for kmers, name in invkmerDict_file1.items():
			if x in kmers:
				for n in name:
					outDict_file1[n]= fastaDict_file1[n]
	outDict_file2= defaultdict(str)
	for x in exclusivetofile2:
		for kmers, name in invkmerDict_file2.items():
			if x in kmers:
				for n in name:
					outDict_file2[n]= fastaDict_file2[n]

	#Checking if "file1_exclusive" and "file2_exclusive" directories exist; if not, create them
	d1= "./file1_exclusive"
	d2= "./file2_exclusive"
	if not os.path.exists(d1) and not os.path.exists(d2):
		os.mkdir(d1)
		os.mkdir(d2)
	
	#Writing out two fasta files to above directories, one with sequences from file 1 that contain kmers not present in file 2 and vice versa
	fp1= d1 + "/uniq_" + os.path.basename(args.file1)	#Given a path, os.path.basename() will extract the file name and return it as a string
	fp2= d2 + "/uniq_" + os.path.basename(args.file2)	#Given a path, os.path.basename() will extract the file name and return it as a string
	ft.write_fasta_dict(outDict_file1, fp1) 
	ft.write_fasta_dict(outDict_file2, fp2)

	#Checking if the combo directories exist; if not, create them
	cd1= "./cmbo_file1_w_uniqfile2"
	cd2= "./cmbo_file2_w_uniqfile1"
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
	cfp1= cd1 + "/cmbo_" + os.path.basename(args.file1)	#Given a path, os.path.basename() will extract the file name and return it as a string
	cfp2= cd2 + "/cmbo_" + os.path.basename(args.file2)	#Given a path, os.path.basename() will extract the file name and return it as a string
	args1= "cp %s %s" % (args.file1, cfp1)
	args2= "cp %s %s" % (args.file2, cfp2)
	subprocess.call(args1, shell=True)
	subprocess.call(args2, shell=True)
	append_fasta_dict(outDict_file2, cfp1)
	append_fasta_dict(outDict_file1, cfp2)
	print("Number of seqs in file 1 [%s]: %d" % (os.path.basename(args.file1), len(fastaDict_file1)))
	print("Number of seqs in file 2 [%s]: %d" % (os.path.basename(args.file2), len(fastaDict_file2)))
	print("Number of seqs with kmers exclusive to file 1: %d" % len(outDict_file1))
	print("Number of seqs with kmers exclusive to file 2: %d" % len(outDict_file2))
	print("Number of seqs in combo fasta v1 (seqs from file 1 + seqs with kmers exclusive to file 2): %d" % (len(fastaDict_file1) + len(outDict_file2)))
	print("Number of seqs in combo fasta v2 (seqs from file 2 + seqs with kmers exclusive to file 1): %d\n" % (len(fastaDict_file2) + len(outDict_file1)))


#If "--subsetMode" flag IS provided
else:
	#Creating an inverse dictionary containing kmers:name
	for key, value in kmerDict_file2.items():
		kmerDict_file2[key]= ", ".join(value)
	invkmerDict_file2 = defaultdict(list)
	{invkmerDict_file2[v].append(k) for k, v in kmerDict_file2.items()}

	#Set of unique kmers for each file, filtering out kmers that contain "X" which are ambiguous amino acids
	uniqkmerSet_file1= kt.kmerSetFasta(args.file1, args.ksize, filter=["X"])
	uniqkmerSet_file2= kt.kmerSetFasta(args.file2, args.ksize, filter=["X"])

	#Any sequences with kmers exclusive to file2 are added to a dictionary
	exclusivetofile2= uniqkmerSet_file2.difference(uniqkmerSet_file1)

	outDict_file2= defaultdict(str)
	for x in exclusivetofile2:
		for kmers, name in invkmerDict_file2.items():
			if x in kmers:
				for n in name:
					outDict_file2[n]= fastaDict_file2[n]

	#Checking if "file2_exclusive" directory exist; if not, create it
	d2= "./file2_exclusive"
	if not os.path.exists(d2):
		os.mkdir(d2)
	
	#Writing out fasta file to above directory, contains sequences from file 2 that contain kmers not present in file 1
	fp2= d2 + "/uniq_" + os.path.basename(args.file2)	#Given a path, os.path.basename() will extract the file name and return it as a string
	ft.write_fasta_dict(outDict_file2, fp2)

	#Checking if the combo fasta directory exists; if not, create it
	cd1= "./cmbo_file1_w_uniqfile2"
	if not os.path.exists(cd1):
		os.mkdir(cd1)

	#Function to append to a fasta file
	def append_fasta_dict(fD, filename):
		fout=open(filename, 'a')
		for k,v in fD.items():
			fout.write(">%s\n%s\n" % (k, v))
		fout.close()

	#Copying file1 to above directory, to which the unique sequences from file2 will be appended to create combo file
	cfp1= cd1 + "/cmbo_" + os.path.basename(args.file1)	#Given a path, os.path.basename() will extract the file name and return it as a string
	args1= "cp %s %s" % (args.file1, cfp1)
	subprocess.call(args1, shell=True)
	append_fasta_dict(outDict_file2, cfp1)
	print("Number of seqs in file 1 [%s]: %d" % (os.path.basename(args.file1), len(fastaDict_file1)))
	print("Number of seqs in file 2 [%s]: %d" % (os.path.basename(args.file2), len(fastaDict_file2)))
	print("Number of seqs with kmers exclusive to file 2: %d" % len(outDict_file2))
	print("Number of seqs in combo fasta (seqs from file 1 + seqs with kmers exclusive to file 2): %d\n" % (len(fastaDict_file1) + len(outDict_file2)))
