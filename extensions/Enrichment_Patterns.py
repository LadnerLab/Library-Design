#!/usr/bin/env python
# coding: utf-8


#import dependancies

import inout as io
import pandas as pd
import glob
from collections import Counter
import argparse


#input arguments and define variables

parser = argparse.ArgumentParser()
parser.add_argument('-d','--input_directory',required=True,help='Input directory with raw enrichment files.')
parser.add_argument('-f','--file_template',required=True,help='Files to be read in (Ex. "XX*.txt".)')
parser.add_argument('-m','--matches',required=True,help='File of protein/peptide matches.')
#parser.add_argument('-l','--lengths',required=True,help='File of protein lengths.')
args=parser.parse_args()

input_directory=args.input_directory
file_template=args.file_template
matches=args.matches
#lengths=args.lengths
#output_file=args.output_file


#read in matches and make dictionary to be "Matches" column.

Matches=io.fileDictHeader(matches, 'CodeName', 'Protein')
Match_prot=Counter(Matches.values())
MatchD=dict(Match_prot)


#read in directory and files; make list with each file.

enrF=glob.glob("%s/%s" % (input_directory,file_template))


#parse enrF to 

counts=[]
peps=[]
enrich=[]
files=[]
for each in enrF:
    files.append(each)
    with open(each,'r')as fin:
        for x in fin:
            for k,v in Matches.items():
                if k in x:
                    pep=x.replace('\n','')
                    peps.append(pep)
                    counts.append(Matches[pep])
    contents=io.fileList(each, col=0, delim='\n', header=False)
    enrich.append(contents)





files_contents=dict(zip(files,(enrich)))





pro_pep={}
count_set=set(counts)
for c in count_set:
    pro_pep[c]=[]





for k,v in pro_pep.items():
    for p in peps:
        if Matches[p]==k:
            v.append(p)





pro_pep_set={}
for pro, pep in pro_pep.items():
    pro_pep_set[pro]=set(pep)





pro_freq={}
for key, value in pro_pep_set.items():
    pro_freq[key]=0
    for file, content in files_contents.items():
        if len(value.intersection(content))!=0:
            pro_freq[key]+=1
        else:
            pro_freq[key]==0





pro_count={}
for pro, pep in pro_pep.items():
    pro_count[pro]=len(pep)

#create dataframe that is exported as txt output file

#Data1=pd.read_csv(lengths, sep='\t',header=0)

Data2=pd.DataFrame(MatchD.items(),columns=["Protein","Matches"])

Data3=pd.DataFrame(pro_count.items(),columns=["Protein","Enrichment"])

Data4=pd.DataFrame(pro_freq.items(),columns=["Protein","Samples"])

#Data12=pd.merge(Data2,Data1)
Data23=pd.merge(Data2,Data3)
Data234=pd.merge(Data23,Data4)

Data234["Enrich/Match"]=Data234["Enrichment"]/Data234["Matches"]
#Data234["Enrich/Length"]=Data234["Enrichment"]/Data234["Length"]
#Data234["Match/Length"]=Data234["Matches"]/Data234["Length"]
Data234["Frequency"]=Data234["Samples"]/len(enrF)

Data234.to_csv(r'output.txt',index=False,header=True,sep='\t')

