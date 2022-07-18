#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import fastatools as ft
import inout as io
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('-u','--unassigned_filename',required=True, help="FASTA file of sequences with unassigned regions.")
parser.add_argument('-m','--metadata',required=True, help="txt file of metadata, see documentation.")
parser.add_argument('-r', '--subject_ratio', type=float, required=True, help='threshold for subject ratio (ex .75)')
args=parser.parse_args()

#Assigns input files for variables.

unassigned_file=args.unassigned_filename
metadata_file=args.metadata
subject_ratio=args.subject_ratio


refD=ft.read_fasta_dict_upper(unassigned_file)
for k in list(refD.keys()):
    if ",," in k:
        refD[k.replace(",,",", ")]=refD[k]






blasts=io.fileDictHeader(metadata_file,"Ref Fasta","Blast Result")    






proteins=io.fileListHeader(metadata_file, "Protein")
proteindict={prot:{}for prot in proteins}






for ref,blast in blasts.items():
    new_fasta=ft.read_fasta_dict_upper(ref)
    with open(blast,'r')as fin:
        fd=pd.read_csv(fin,sep='\t',header=0)
        fd["Subject_Ratio"]=(fd["Subject End"]-fd["Subject Start"]+1)/fd["Subject Length"]
        fd=fd[fd["Subject_Ratio"]>=.67]
        for p in proteins:
            if p in blast:
                fd["Protein"]=p
        fd["Query&Protein"]=fd["Query Name"]+ " "+fd["Protein"]
        fd["Query&Range"]=fd["Query Name"]+ " "+fd["Query Start"].apply(str)+'-'+fd["Query End"].apply(str)
        fd2=fd[["Query&Protein","Query Name","Query Start","Query End"]]
        fd3=fd[["Query&Range","Query Name","Query Start","Query End"]]

        combo=fd2.values.tolist()
        combo2=fd3.values.tolist()
        
        names=[]
        seqs=[]
        for x in combo2:
            names.append(x[0])
            seqs.append(refD[x[1]][x[2]:x[3]])
            
        new_entry=dict(zip(names,seqs))
       
        for p in proteins:
            for index, row in fd.iterrows():
                if p in row["Protein"]:
                    for names, seqs in new_entry.items():
                        proteindict[p][names]=[seqs]
                        new_fasta.update(proteindict[p])
                        ft.write_fasta_dict(proteindict[p],'%s.fasta'%p)
            
                        







