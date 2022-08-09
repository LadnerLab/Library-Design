#!/usr/bin/env python3

import argparse, os	
from collections import defaultdict

#Example command: cluster_dir_sort.py -c . -f 40119_Parvovirinae_clusterinfo_Complete.tsv -p 40119_id_70_ 

parser = argparse.ArgumentParser(description='''After manually assigning each cluster to a protein within the second column of the output file from 
					cluster_annotations_parse.py, this script will create directories for each protein and move the corresponding clusters to 
					their appropriate protein directory based off of the assignments.''')

parser.add_argument("-p", "--prefix", default="", metavar='\b', help="File name prefix of clusters preceding cluster number, if applicable (e.g. 11158_id_70_).")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-f", "--file", required=True, metavar='\b', help='Name of output file from cluster_annotations_parse.py, containing cluster numbers and manually added protein assignments.')
reqArgs.add_argument("-c", "--clusters", required=True, metavar='\b', help="Path to directory containing clusters for sorting. Do not include slash at end of path. Note, cluster names must end with cluster number.")

args = parser.parse_args()


with open(args.file, "r") as fin:
	protLabel={}
	next(fin)		#Skips header line
	for line in fin:
		l= line.rstrip("\n").split("\t")
		protLabel[l[0]]=l[1]

#Creating directory for each unique protein assignment
for p in set(protLabel.values()):
	dirP= args.clusters + "/" + p
	if not os.path.exists(dirP):
		os.mkdir(dirP)

for clusterNum, prot in protLabel.items():
	if len(prot) == 0:
		print("Cluster %s has not been assigned to a protein. Cannot be sorted into a subdirectory." % clusterNum)
	
	#Moving clusters to corresponding protein directory based off of assignments in input file (column 2)
	fileName = args.prefix + clusterNum
	filePath= args.clusters + "/" + fileName
	fileDest= args.clusters + "/" + prot + "/" + fileName
	os.rename(filePath, fileDest)