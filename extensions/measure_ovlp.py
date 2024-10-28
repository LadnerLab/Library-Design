#!/usr/bin/env python3
import os
import pandas as pd
import kmertools as kt
import argparse
import concurrent.futures

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i1','--HV2-metadata', help='Filepath to tab-delimited HV2 metadata', required=True)
    parser.add_argument('-i2','--HV3-metadata', help='Filepath to tab-delimited HV3 metadata', required=True)
    parser.add_argument('-ks', '--kmer-size', default=6, help='Filepath to tab-delimited HV3 metadata', required=False)
    parser.add_argument('-m', '--multiprocessed', action="store_true", help='Run analysis with multiple processes', required=False)
    parser.add_argument('-o', '--output-dir', help='Output directory to output peptides and overlap values', required=True)
    
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # extract data
    HV2_df = pd.read_csv(args.HV2_metadata, usecols=["CodeName", "Species", "Peptide"], sep="\t")
    HV3_df = pd.read_csv(args.HV3_metadata, usecols=["CodeName", "Species", "Peptide"], sep="\t")

    # if HV3 Species column is split with ';', then only use first one
    HV3_df["Species"] = HV3_df["Species"].str.split(';').str[0]

    # get max overlap for each peptide
    overlap_data = get_overlap_species_by_species(HV2_df, HV3_df, args.kmer_size, args.multiprocessed)

    # save the dataframe
    out_df = pd.DataFrame(overlap_data, columns=["HV2_CodeName", "HV3_CodeName", "MaxOvlpScore"])
    out_df.to_csv(os.path.join(args.output_dir, "peptide_ovlp_scores.tsv"), sep='\t', index=False)


def get_overlap_species_by_species(HV2_df:pd.DataFrame, HV3_df:pd.DataFrame, kmer_size:int, multiprocessed:bool)->list:
    # [(HV2_CodeName, HV3_CodeName, MaxOvlpScore), (HV2_CodeName, HV3_CodeName, MaxOvlpScore), ...]
    overlap_data = list()
    HV2_species_to_HV3_species = dict()

    # map the species names between the two
    HV2_species_set = set(HV2_df["Species"].tolist())
    HV3_species_set = set(HV3_df["Species"].tolist())
    for HV2_species in HV2_species_set:
        # get the corresponding HV3 species
        mapped_species_list = [HV3_species for HV3_species in HV3_species_set if HV2_species.replace(' ', '_') in HV3_species]

        # map the species if one was found
        if len(mapped_species_list) > 0:
            HV2_species_to_HV3_species[HV2_species] = mapped_species_list[0]
        else:
            # let user know that none was found for species
            print(f"No HV3 species found for {HV2_species}")

    # get overlap for each species
    HV2_species_groups = HV2_df.groupby("Species")
    HV3_species_groups = HV3_df.groupby("Species")

    if multiprocessed:
        # run each species in a different process
        with concurrent.futures.ProcessPoolExecutor() as process_executor:
            futures = [process_executor.submit( get_overlap_single_species,
                                            HV2_species_groups.get_group(HV2_species), 
                                            HV3_species_groups.get_group(HV3_species),
                                            kmer_size,
                                            HV2_species
                                            ) for HV2_species, HV3_species in HV2_species_to_HV3_species.items()]

            for future in concurrent.futures.as_completed(futures):
                single_species_overlap_list = future.result()

                # concat list
                overlap_data += single_species_overlap_list
    else:
        for HV2_species, HV3_species in HV2_species_to_HV3_species.items():
            overlap_data += get_overlap_single_species(HV2_species_groups.get_group(HV2_species), HV3_species_groups.get_group(HV3_species), kmer_size, HV2_species)
             
    return overlap_data


def get_overlap_single_species(HV2_df:pd.DataFrame, HV3_df:pd.DataFrame, kmer_size:int, species:str)->list:
    print(f"Working on {species}")

    overlap_data = list()

    # create kmer_sets
    HV2_kmer_sets = dict()
    HV3_kmer_sets = dict()

    for i, HV2_row in HV2_df.iterrows():
        HV2_codename = HV2_row["CodeName"]
        HV2_peptide = HV2_row["Peptide"]
        HV2_kmer_sets[HV2_codename] = kt.kmerSet(HV2_peptide,kmer_size)

    for i, HV3_row in HV3_df.iterrows():
        HV3_codename = HV3_row["CodeName"]
        HV3_peptide = HV3_row["Peptide"]
        HV3_kmer_sets[HV3_codename] = kt.kmerSet(HV3_peptide,kmer_size)

    # get max overlap for each peptide
    for HV2_codename, HV2_kmer_set in HV2_kmer_sets.items():
        max_ovlp_score = 0
        max_HV3_codename = str()

        for HV3_codename, HV3_kmer_set in HV3_kmer_sets.items():
            # get the overlap score between the peptides
            ovlp_score = get_overlap(HV2_kmer_set, HV3_kmer_set)

            # set it to max for this HV2 peptide if it is the greatest
            if ovlp_score >= max_ovlp_score:
                max_ovlp_score = ovlp_score
                max_HV3_codename = HV3_codename

        overlap_data.append((HV2_codename, max_HV3_codename, max_ovlp_score))

    return overlap_data


def get_overlap(HV2_kmer_set:set, HV3_kmer_set:set)->float:
    # get the overlap (percentage of kmers that overlap)
    ovlp_score = len(HV2_kmer_set.intersection(HV3_kmer_set))/len(HV2_kmer_set)

    return ovlp_score


if __name__ == "__main__":
    main()

