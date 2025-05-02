#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
import fastatools as ft

# This scripts generates a dataset to use as training data for the oligo_encoding neural network that predicts bindings

NTS = ['A', 'C', 'G', 'T']
AAS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
CODONS = [f'{n1}{n2}{n3}' for n1 in NTS for n2 in NTS for n3 in NTS]
#CODONS = ['AAA', 'AAC', 'ACC', 'ACT', 'ATC', 'CAC', 'CAG', 'CCG', 'CGC', 'CGT', 'CTG', 'GAA', 'GAC', 'GCT', 'GGC', 'GGT', 'GTT', 'TAC', 'TCC', 'TCT', 'TGC', 'TTC']
'''
AA_CATEGORIES = {
    "aromatic": {'F', 'W', 'Y'},
    "nonpolar": {'A', 'G', 'I', 'L', 'M', 'V'},
    "polar": {'P', 'C', 'N', 'Q', 'S', 'T'},
    "positive charged": {'H', 'K', 'R'},
    "negative charged": {'D', 'E'}
}
'''

# print(f"Num Columns: {len(NTS)+len(AAS)+len(CODONS)+len(AA_CATEGORIES)}")
print(f"Num Columns: {len(NTS)+len(AAS)+len(CODONS)}")

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('-m', '--metadata', help='Tab-delimited file with columns for ' \
            'Library Name, Reads Path (if training set), NN Path, and AA Path', required=True)
    reqArgs.add_argument('-o', '--output_file', help='Filename to output training matrix.', required=True)

    optArgs = parser.add_argument_group('Optional Arguments')
    optArgs.add_argument('-z', '--include_score', action='store_true', help='Include this file if Reads Path is ' \
        'included in metadata (i.e., for generating a training matrix).', required=False)
    optArgs.add_argument('-a', '--average_samples', action='store_true', help='Average counts between ' \
            'replicate samples with the same basename (everything before the last underscore) before calculating zscores.', required=False)

    args = parser.parse_args()

    # create/reset the output file
    reset_file(args.output_file, args.include_score)

    # parse through each library in metadata and generate matrix
    metadata_df = pd.read_csv(args.metadata, sep="\t", index_col=0)
    for i, row in metadata_df.iterrows():
        counts_file = None
        if args.include_score:
            counts_file = row["Reads Path"]
        nt_file = row["NN Path"]
        aa_file = row["AA Path"]

        matrix_df = gen_feature_matrix_single_lib(counts_file, nt_file, aa_file, args.average_samples)
        
        # append to file
        matrix_df.to_csv(args.output_file, mode='a', sep="\t", header=False, index=False)


def gen_feature_matrix_single_lib(counts_file, nt_file, aa_file, average_samples):
    # read in data
    aa_seqs = ft.read_fasta_dict(aa_file)
    nt_seqs = ft.read_fasta_dict(nt_file)

    # generate ratios (indexed on sequence name)
    ratio_df = generate_ratios(aa_seqs, nt_seqs)

    if counts_file:
        counts_df = pd.read_csv(counts_file, sep="\t", index_col=0)

        if average_samples:
            # average counts for sampels with the same basename
            counts_df = average_sample_counts(counts_df)
        
        # convert raw_counts to scores
        scores_df = calculate_zscores(counts_df)

        # ensure that most dataframes have the same index values
        ratio_df, scores_df = interect_index(ratio_df, scores_df)

        # generate neural netword trainging matrix
        nn_df = generate_training_matrix(scores_df, ratio_df)
    else:
        nn_df = ratio_df.reset_index(drop=False)
    
    return nn_df


def generate_training_matrix(scores_df, ratio_df):
    out_data = list()

    for seq_name, row in ratio_df.iterrows():
        row = row.to_list()
        row.insert(0, np.nan)
        row.insert(1, np.nan)
        for sample_name in scores_df.columns.to_list():
            temp_row = row.copy()
            score = scores_df.loc[seq_name, sample_name]
            temp_row[0] = seq_name
            temp_row[1] = score
            out_data.append(temp_row)
    
    out_df = pd.DataFrame(out_data, columns=["Sequence name", "Score"]+ratio_df.columns.to_list())

    return out_df

def generate_ratios(aa_seqs, nt_seqs):
    encoding_data = dict()
    for name in nt_seqs.keys():
        encoding_data[name] = {
            "nt_counts": {nt: 0 for nt in NTS},
            "aa_counts": {aa: 0 for aa in AAS},
            "codon_counts": {codon: 0 for codon in CODONS}
        }

    found_peps = set()
    found_aa_data = dict()
    
    for name, seq in nt_seqs.items():
        # nt counts and codon counts
        for i in range(len(seq)):
            nt = seq[i]
            encoding_data[name]["nt_counts"][nt] += 1

            if i < len(seq) - 2:
                codon = seq[i:i+3]
                if codon in CODONS:
                    encoding_data[name]["codon_counts"][codon] += 1
    
        # aa counts
        pep_basename = name.split("-")[0]
        if pep_basename in found_peps:
            encoding_data[name]["aa_counts"] = found_aa_data[pep_basename]
        else:
            for aa in aa_seqs[pep_basename]:
                encoding_data[name]["aa_counts"][aa] += 1
            found_peps.add(pep_basename)
            found_aa_data[pep_basename] = encoding_data[name]["aa_counts"]

    ratio_df = create_ratio_df_from_dict(encoding_data)

    # complexity
    # ratio_df["Complexity"] = len(nt_seqs)
    print(ratio_df)

    return ratio_df

def create_ratio_df_from_dict(encoding_data):
    df_data = list()

    for name, counts_dict in encoding_data.items():
        row_data = [name]

        # nt ratio
        row_data += get_ratio_data(list(counts_dict["nt_counts"].values()))

        # aa ratio
        row_data += get_ratio_data(list(counts_dict["aa_counts"].values()))

        # codon ratio
        row_data += get_ratio_data(list(counts_dict["codon_counts"].values()))

        # categorical ratio
        ''''
        aa_seq_len = sum(list(counts_dict["aa_counts"].values()))
        for aa_set in AA_CATEGORIES.values():
            row_data.append(sum([counts_dict["aa_counts"][aa] for aa in aa_set]) / aa_seq_len)
        '''

        df_data.append(row_data)
    
    # return pd.DataFrame(df_data, columns=["Sequence name"]+NTS+AAS+CODONS+list(AA_CATEGORIES.keys())).set_index("Sequence name")
    return pd.DataFrame(df_data, columns=["Sequence name"]+NTS+AAS+CODONS).set_index("Sequence name")

def get_ratio_data(counts_list):
    ret_data = list()
    total = sum(counts_list)
    for count in counts_list:
        ret_data.append(count/total)
    return ret_data

def interect_index(df1, df2):
    common_index = df1.index.intersection(df2.index)
    return df1.loc[common_index], df2.loc[common_index]

def calculate_zscores(counts_df):
    log_transform_df=np.log(counts_df+1)
    zscore_transform_df=(log_transform_df-log_transform_df.mean())/log_transform_df.std()
    return zscore_transform_df

def reset_file(output_file, include_reads):
    with open(output_file, "w") as f:
        ''' No column names
        columns = NTS+AAS+CODONS
        if include_reads:
            columns = ["Sequence name", "Sample name", "Score"]+NTS+AAS+CODONS
        columns_str = '\t'.join(columns)
        f.write(columns_str+"\n")
        '''
        pass

def average_sample_counts(counts_df):
    average_counts_df = pd.DataFrame(index=counts_df.index)
    replicate_groups = get_sample_replicates(counts_df.columns.to_list())
    for base, replicates in replicate_groups.items():
        average_counts_df[base] = counts_df[replicates].mean(axis=1)
    return average_counts_df

def get_sample_replicates(samples):
    seen_basenames = set()
    sample_pairs = dict()

    for sample in samples:
        basename = get_basename(sample)
        if basename in seen_basenames:
            sample_pairs[basename].append(sample)
        else:
            seen_basenames.add(basename)
            sample_pairs[basename] = [sample]

    return sample_pairs

def get_basename(name):
    return "_".join(name.split("_")[:-1])

if __name__ == "__main__":
    main()