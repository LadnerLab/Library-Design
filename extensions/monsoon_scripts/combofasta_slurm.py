#!/usr/bin/env python3

import argparse, os
from collections import defaultdict
from subprocess import Popen, PIPE


parser = argparse.ArgumentParser(description='''A wrapper script for batch processing pairs of files listed in a metafile with combofasta.py on slurm.''')

parser.add_argument("-k", "--ksize", default=9, type=int, metavar='\b', help="kmer size to use for comparison [default: 9]")
parser.add_argument("-s", "--subsetMode", default=False, action="store_true", help="This flag can be used if you only wish to check if file2 is a subset of file1. In this mode, only one combo fasta file version will be created— seqs from file1 + seqs w kmers unique to file2. A seperate fasta file containing only the latter will be created as well (file2 seqs that contain kmers not present in file1).")
parser.add_argument("-p", "--path", default="combofasta.py", metavar='\b', help="Path of combofasta.py [default: combofasta.py]")
#Resource utilization
parser.add_argument("-c", "--cpus", default=1, type=int, metavar='\b', help="Use this to specify the number of cpus per task [default: 1].")
parser.add_argument("-m", "--mem", default=1, type=int, metavar='\b', help="Use this to specify the amount of memory to allow be used, in G [default: 1].")
parser.add_argument("-t", "--time", default="20:00", metavar='\b', help="Use this to specify the amount of time needed for script to run [default: 20:00].")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-mf", "--metafile", required=True, metavar='\b', help="Tab-delimited file with 2 columns, with each row containing a pair of files to be compared.")
reqArgs.add_argument("-f1", "--file1", choices=[0,1], type=int, required=True, metavar='\b', help="Column number of fasta file #1 in metafile [choices: 0,1]")
reqArgs.add_argument("-f2", "--file2", choices=[0,1], type=int, required=True, metavar='\b', help="Column number of fasta file #2 in metafile [choices: 0,1]")

args = parser.parse_args()


#Reading in metafile and storing names of file pairs in a dictionary
metaDict = defaultdict(str)
with open(args.metafile, "r") as fin:
	for line in fin:
		line= line.strip().split("\t")
		metaDict[line[args.file1]]= line[args.file2]

#Checking if directories for stdout and stderr exist; if not, create them
do= "./stdout"
de= "./stderr"
if not os.path.exists(do) and not os.path.exists(de):
	os.mkdir(do)
	os.mkdir(de)

for file1, file2 in metaDict.items():
	if args.subsetMode:
		cmd= "%s -f1 %s -f2 %s -k %d --subsetMode" % (args.path, file1, file2, args.ksize)
	else:
		cmd= "%s -f1 %s -f2 %s -k %d" % (args.path, file1, file2, args.ksize)
	f1name=	os.path.basename(file1)	#Given a path, os.path.basename() will extract the file name and return it as a string
	f2name= os.path.basename(file2)	#Given a path, os.path.basename() will extract the file name and return it as a string
	n= os.path.splitext(f1name)[0] + "-" + os.path.splitext(f2name)[0]
	fpo= do + "/" + n + ".out"
	fpe= de + "/" + n + ".err"
	cmd= "sbatch -J cmbofaa-%s -o %s -e %s -t %s -c %d --mem=%dG --wrap='%s'" % (n, fpo, fpe, args.time, args.cpus, args.mem, cmd)
	print (cmd)
	Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
	