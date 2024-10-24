#!/usr/bin/env python3
import glob
import os
import kmertools as kt
import fastatools as ft
import itertools as it
import inout as io
from collections import defaultdict
import argparse
import subprocess
import pandas as pd


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i', '--batch_map', required=True,
            help="Tab-delimited file mapping target alignments file to subtype directory.")
    parser.add_argument('-k', '--kmer-size', type=int, default=6, help='', required=False)
    parser.add_argument('--kmer-ovlp-thresh', type=float, default=0.3, help='')
    parser.add_argument('-o', '--output-dir',  help='Output directory for where to output protein combination files.', required=True)
    
    args = parser.parse_args()

    assert not os.path.exists(args.output_dir), \
        f"Warning: The directory \"{args.output_dir}\" already exists. Please move or delete."

    os.mkdir(args.output_dir)

    # create map from target alignment file to subtype directory (checks for valid filepaths)
    batch_map = create_batch_map(args.batch_map)

    generateSubtypeAlignmentsMultipleSpecies(
            batch_map=batch_map, 
            output_dir=args.output_dir, 
            ks=args.kmer_size, 
            kmer_ovlp_thresh=args.kmer_ovlp_thresh
            )

def generateSubtypeAlignmentsMultipleSpecies(batch_map, output_dir, ks, kmer_ovlp_thresh):
    # create protein assignments/alignments for each batch
    for target_alignments_file, subtype_dir in batch_map.items():
        subtype_out_dir = make_dir(output_dir, f"{os.path.basename(subtype_dir)}_ByProtein")

        generateSubtypeAlignmentsSingleSpecies(
                target_file=target_alignments_file, 
                subtype_dir=subtype_dir, 
                ks=ks, 
                kmer_ovlp_thresh=kmer_ovlp_thresh, 
                output_dir=subtype_out_dir
                )


def create_batch_map(batch_map_filepath):
    batch_map = dict(pd.read_csv(batch_map_filepath, sep='\t').values)

    # make sure all entries exist
    for target_alignments_file, subtype_dir in batch_map.items():
        if not os.path.exists(target_alignments_file):
            raise FileNotFoundError(f"{target_alignments_file} does not exist.")
        elif not os.path.exists(subtype_dir):
            raise FileNotFoundError(f"{subtype_dir} does not exist.")

    return batch_map


def generateSubtypeAlignmentsSingleSpecies(target_file, subtype_dir, ks, kmer_ovlp_thresh, output_dir):
    # dataframe for mapping protein to file
    target_files_data = list()

    # create alignment dir
    align_dir = make_dir(output_dir, "aligned")

    targetD = io.fileDictHeader(target_file, "Protein", "Fasta")
    targetKmersD = {k:kt.kmerSetFasta(v, ks) for k,v in targetD.items()}

    proteinD = {pn:defaultdict(dict) for pn in targetD}

    fastaL = glob.glob(os.path.join(subtype_dir, "*fasta"))

    # assign seqs to proteins
    for each in fastaL:
        fD = ft.read_fasta_dict(each)
        fD = rmvRedundant(fD)

        # keep track of overlap scores to only add the sequence with the greatest overlap
        ovlp_dict = {pn:kmer_ovlp_thresh for pn in targetD}
        seqs_to_include = {pn:None for pn in targetD}

        for n,s in fD.items():
            kmers = kt.kmerSet(s,ks)
            kOvlp = {p:len(kmers.intersection(pk))/len(kmers) for p,pk in targetKmersD.items()}
            if max(kOvlp.values())>=kmer_ovlp_thresh:
                topProt = sorted([(v,k) for k,v in kOvlp.items()])[::-1][0][1]

                # check if this is highest overlap score for this protein
                if max(kOvlp.values()) >= ovlp_dict[topProt]:
                    # assign new sequence to use from this fasta
                    seqs_to_include[topProt] = (n, s)
                    ovlp_dict[topProt] = max(kOvlp.values())

            '''
            else:
                print(kOvlp)
            '''

        # assign sequences to proteins
        for topProt, seq_data in seqs_to_include.items():
            # test if a sequence was assigned for this protein
            if seq_data != None:
                proteinD[topProt][seq_data[0]] = seq_data[1]

    for pn, fD in proteinD.items():
        protein_fasta_path = os.path.join(output_dir, f"{pn}_combo.fasta")
        ft.write_fasta_dict(fD, protein_fasta_path)

        # warn user if multiple protein sequences were assigned to one cluster
        if len(fD) < 1:
            print(f"Warning: no protein sequences from {os.path.basename(subtype_dir)} were assigned to {pn}")

        align_file(
            in_fasta=protein_fasta_path,
            out_fasta=os.path.join(align_dir, f"einsi_{pn}_combo.fasta"),
            quiet=True
            )

        # append protein and fasta
        target_files_data.append((pn, os.path.join(align_dir, f"einsi_{pn}_combo.fasta")))

    # save target files fasta
    target_file = pd.DataFrame(target_files_data, columns=["Protein", "Fasta"])
    target_file.to_csv(os.path.join(output_dir, f"target_files.tsv"), sep='\t', index=False)


def isSubString(s1, s2):
    longerKs = kt.kmerSet(s1,len(s2))
    if s2 in longerKs:
        return True
    else:
        return False


# Remove sequences that are entirely redundant of other sequences in the same input dictionary
def rmvRedundant(fD, verbose=False):
    toExclude=defaultdict(list)
    seqByLen = [s[1] for s in sorted([(len(v), k) for k,v in fD.items()])[::-1]]
    for a,b in it.combinations(seqByLen, 2):
        if isSubString(fD[a],fD[b]):
            toExclude[b].append(a)
    outD = {n:s for n,s in fD.items() if n not in toExclude}
    if verbose and len(toExclude)>0:
        print(toExclude)
    return outD


# align file with mafft
def align_file(in_fasta, out_fasta, quiet):
    command = ["mafft-einsi", "--preservecase", "--inputorder", "--thread", "-1", in_fasta]
    if quiet:
        command.insert(1, "--quiet")

    result = subprocess.run(command, stdout = subprocess.PIPE, universal_newlines = True)

    with open(out_fasta, 'w') as out_fasta:
        # remove newline characters that mafft adds
        lines = result.stdout.split('\n')
        for line_idx in range(len(lines)):
            line = lines[line_idx]
            if line_idx == 0:
                out_fasta.write(line + '\n')
            elif line.startswith('>'):
                out_fasta.write('\n' + line + '\n')
            else:
                out_fasta.write(line)


def make_dir(path, new):
    dir_name = os.path.join(path, new)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    return dir_name


if __name__ == "__main__":
    main()