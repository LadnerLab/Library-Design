#!/usr/bin/env python3
import os
import glob
import pandas as pd
import fastatools as ft
import argparse

from BatchAssign import generateSubtypeAlignmentsSingleSpecies
from CreateBridge import create_bridge_alignments_single_species

WINDOW_SIZE = 30


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i','--file-metadata', help='', required=True)
    parser.add_argument('-k', '--kmer-size', type=int, default=6, help='Used for assigning sequences to proteins. Kmer size for kmer sets that are going to be intersected.', required=False)
    parser.add_argument('--kmer-ovlp-thresh', type=float, default=0.3, help='Used for assigning sequences to proteins. Minimum kmer overlap that the largest overlapping sequence must have to be assigned to a protein.', required=False)
    parser.add_argument('--max-window-gaps', type=int, default=WINDOW_SIZE/2, required=False)

    parser.add_argument('-o', '--output-dir', 
        help='Output directory for where to output protein combination files and bridge alignments.', required=True)
    
    args = parser.parse_args()

    assert not os.path.exists(args.output_dir), \
        f"Warning: The directory \"{args.output_dir}\" already exists. Please move or delete."

    os.mkdir(args.output_dir)

    map_bridge_and_link(
        file_metadata=args.file_metadata,
        kmer_size=args.kmer_size,
        kmer_ovlp_thresh=args.kmer_ovlp_thresh,
        max_window_gaps=args.max_window_gaps,
        output_dir=args.output_dir
        )


def map_bridge_and_link(
    file_metadata,
    kmer_size,
    kmer_ovlp_thresh,
    max_window_gaps,
    output_dir
    ):
    # create map for assigning seqs to proteins in batches
    batch_map = create_filepath_map(file_metadata)

    # create protein assignments/alignments for each batch
    for i, row in batch_map.iterrows():
        # extract data from columns
        output_name = row["Name"]
        target_alignments_file = row["TargetMapFile"]
        subtype_dir = row["SubtypeDir"]
        bridge_alignments_file = row["BridgeMapFile"]
        epitope_positions_file = row["EpitopeFile"]

        map_bridge_and_link_single_row(
            kmer_size=kmer_size,
            kmer_ovlp_thresh=kmer_ovlp_thresh,
            max_window_gaps=max_window_gaps,
            output_dir=output_dir,
            output_name=output_name,
            target_alignments_file=target_alignments_file,
            subtype_dir=subtype_dir,
            bridge_alignments_file=bridge_alignments_file,
            epitope_positions_file=epitope_positions_file
            )


def map_bridge_and_link_single_row(
    kmer_size,
    kmer_ovlp_thresh,
    max_window_gaps,
    output_dir,
    output_name,
    target_alignments_file,
    subtype_dir,
    bridge_alignments_file,
    epitope_positions_file
    ):

    spec_output_dir = make_dir(output_dir, output_name)

    subtype_out_dir = make_dir(spec_output_dir, f"{os.path.basename(subtype_dir)}_ByProtein")

    # create new species by mapping sequences to proteins from old species
    generateSubtypeAlignmentsSingleSpecies(
            target_file=target_alignments_file, 
            subtype_dir=subtype_dir, 
            ks=kmer_size, 
            kmer_ovlp_thresh=kmer_ovlp_thresh, 
            output_dir=subtype_out_dir
            )

    # TODO: let user defile target files output
    # assume target files for new species
    new_species_target_file = os.path.join(subtype_out_dir, "target_files.tsv")

    mapping_name = f"{bridge_alignments_file.split(os.sep)[-2]}_{new_species_target_file.split(os.sep)[-2]}"
    bridge_out_dir = make_dir(spec_output_dir, f"Bridge_{output_name}_{os.path.basename(subtype_dir)}_ByProtein")

    # create bridge fasta files with a representative from each species
    create_bridge_alignments_single_species(
            old_species_file = bridge_alignments_file,
            new_species_file = new_species_target_file, 
            output_dir = bridge_out_dir
            )

    # assume target files for bridge
    bridge_target_files = os.path.join(bridge_out_dir, "target_files.tsv")

    # read epitope positions
    epitope_positions_df = pd.read_csv(epitope_positions_file, sep='\t')[["ClusterID", "PeptideID", "Start Position", "Stop Position"]]

    out_df = pd.DataFrame(columns=["ClusterID", "PeptideID", "SequenceName", "Window", "Start Position", "Stop Position"])
    removed_df = pd.DataFrame(columns=["ClusterID", "PeptideID", "SequenceName", "Window", "Start Position", "Stop Position"])

    # extract align maps
    # r2a stands for raw to aligned and a2r stands for aligned to raw
    # IMPORTANT: assume that epitope positions were taken from targets align map
    targets_align_map_a2r = extract_align_map( bridge_alignments_file, rev=True )

    bridge_align_map_r2a = extract_align_map( bridge_target_files )
    bridge_align_map_a2r = extract_align_map( bridge_target_files, rev=True )

    subtype_align_map_r2a = extract_align_map( new_species_target_file )
    subtype_fasta_dict = extract_align_map(new_species_target_file, read_fasta=True)

    '''
    # make sure proteins are the same for smalllest
    assert targets_align_map_a2r.keys() == bridge_align_map_r2a.keys() and bridge_align_map_r2a.keys() == subtype_align_map_r2a.keys(), \
        "Protein names do not match for all alignments."
    '''
    # use smallest protein names
    proteins = bridge_align_map_r2a.keys()

    assert all( protein in set(epitope_positions_df["ClusterID"].to_list()) for protein in proteins), \
        "Alignments include a protein name that is not defined in the epitope positions file."

    out_row_idx = 0
    # loop through each protein
    for protein in proteins:
        out_df, removed_df = link_alignments_single_protein(
                                    targets_align_map_a2r, 
                                    bridge_align_map_r2a, 
                                    bridge_align_map_a2r, 
                                    subtype_align_map_r2a, 
                                    out_df, 
                                    protein,
                                    epitope_positions_df,
                                    subtype_fasta_dict[protein],
                                    max_window_gaps,
                                    removed_df
                                    )

    out_df.to_csv(os.path.join(spec_output_dir, f"new_{os.path.basename(epitope_positions_file)}"), sep='\t', index=False)
    removed_df.to_csv(os.path.join(spec_output_dir, f"removed_peptides.tsv"), sep='\t', index=False)


def link_alignments_single_protein(
        targets_align_map_a2r, 
        bridge_align_map_r2a, 
        bridge_align_map_a2r, 
        subtype_align_map_r2a, 
        out_df, 
        protein, 
        epitope_positions_df, 
        fasta_dict,
        max_window_gaps,
        removed_df
        ):

    # get bridge candidates (target candidate is first and subtype candidate is second)
    bridge_target_candidate_name = list(bridge_align_map_r2a[protein].keys())[0]
    bridge_subtype_candidate_name = list(bridge_align_map_r2a[protein].keys())[1]

    target_candidate_name = bridge_target_candidate_name
    subtype_candidate_name = bridge_subtype_candidate_name

    # test if the target name is in the range to be truncated in mafft alignment
    if len(target_candidate_name) > 225:
        # test if target candidate name was truncated in alignment
        if target_candidate_name not in targets_align_map_a2r[protein].keys():
            # replace it with the verse name in the aligned fasta files
            for name in targets_align_map_a2r[protein].keys():
                if target_candidate_name == name[0:len(target_candidate_name)]:
                    target_candidate_name = name

    # test if the subtype name is in the range to be truncated in mafft alignment
    if len(subtype_candidate_name) > 225:
        # test if subtype candidate name was truncated in alignment
        if subtype_candidate_name not in subtype_align_map_r2a[protein].keys():
            # replace it with the verse name in the aligned fasta files
            for name in subtype_align_map_r2a[protein].keys():
                if subtype_candidate_name == name[0:len(subtype_candidate_name)]:
                    subtype_candidate_name = name

    # extract dict data
    target_a2r = targets_align_map_a2r[protein][target_candidate_name]
    bridge_r2a = bridge_align_map_r2a[protein][bridge_target_candidate_name]
    bridge_a2r = bridge_align_map_a2r[protein][bridge_subtype_candidate_name]
    subtype_r2a = subtype_align_map_r2a[protein][subtype_candidate_name]

    # cut to df to only have protein
    cut_df = epitope_positions_df[epitope_positions_df["ClusterID"] == protein]

    peptide_positions_data = list()
    removed_peptides = list()

    # loop through each peptide
    for i, row in cut_df.iterrows():
        target_start_pos = row["Start Position"]
        target_stop_pos = row["Stop Position"]

        # convert start and stop positions

        # convert target aligned to raw
        subtype_start_pos = convert_target_2_subtype(target_start_pos, target_a2r, bridge_r2a, bridge_a2r, subtype_r2a, protein)
        subtype_stop_pos = convert_target_2_subtype(target_stop_pos, target_a2r, bridge_r2a, bridge_a2r, subtype_r2a, protein)

        # print(f"old: ({target_start_pos}, {target_stop_pos})")
        # print(f"new: ({subtype_start_pos}, {subtype_stop_pos})")

        # shift positions for each sequence
        for seq_name, seq in fasta_dict.items():
            # create copies of start and stop for sequence
            seq_start_pos = subtype_start_pos
            seq_stop_pos = subtype_stop_pos


            seq_positions = subtype_align_map_r2a[protein][seq_name]
            # get end position in sequence
            max_pos = list(seq_positions.values())[len(seq_positions) - 1]

            # extract window sequence
            window_seq = seq[seq_start_pos-1:seq_stop_pos-1]

            # check if the window is not entirely out of the sequence, and still a valid peptide
            if seq_stop_pos-seq_start_pos > 1:
                # shift until the number of characters that are not - is length 30
                while( window_seq.count('-') < max_window_gaps and len(window_seq) - window_seq.count('-') < 30 and ( seq_start_pos > 1 or seq_stop_pos < max_pos ) ):
                    # test if shift left is possible
                    if( seq_start_pos - 1 > 0 ):
                        # shift left
                        seq_start_pos -= 1

                    # extract new window sequence
                    window_seq = seq[seq_start_pos-1:seq_stop_pos-1]

                    # test until the number of characters that are not - is length 30 and shift right is possible
                    if( seq_stop_pos-seq_start_pos != 0 and len(window_seq) - window_seq.count('-') != 30 and seq_stop_pos + 1 <= max_pos ):
                        # shift right
                        seq_stop_pos += 1

                        # extract new window sequence
                        window_seq = seq[seq_start_pos-1:seq_stop_pos-1]

                # shift inwards until number of characters is exactly 30
                while len(window_seq) - window_seq.count('-') > 30:
                    # shift right
                    seq_stop_pos -= 1

                    # extract new window sequence
                    window_seq = seq[seq_start_pos-1:seq_stop_pos-1]

                    if len(window_seq) - window_seq.count('-') > 30:
                        # shift left
                        seq_start_pos += 1

                        # extract new window sequence
                        window_seq = seq[seq_start_pos-1:seq_stop_pos-1]

            # add data to df
            # only add if a comparable peptide was chose (check number of gaps)
            if( window_seq.count('-') < max_window_gaps and seq_stop_pos-seq_start_pos > 1):
                peptide_positions_data.append([protein, row["PeptideID"], seq_name, window_seq, seq_start_pos, seq_stop_pos, target_candidate_name, subtype_candidate_name])
            else:
                removed_peptides.append([protein, row["PeptideID"], seq_name, window_seq, seq_start_pos, seq_stop_pos, target_candidate_name, subtype_candidate_name])


    # create df for protein and append to out df
    tmp_df = pd.DataFrame(removed_peptides, columns=["ClusterID", "PeptideID", "SequenceName", "Window", "Start Position", "Stop Position", "TargetCandidate", "SubtypeCandidate"])
    removed_df = pd.concat([removed_df,tmp_df])
    tmp_df = pd.DataFrame(peptide_positions_data, columns=["ClusterID", "PeptideID", "SequenceName", "Window", "Start Position", "Stop Position", "TargetCandidate", "SubtypeCandidate"])
    out_df = pd.concat([out_df,tmp_df])

    # return updated table and current row idx
    return out_df, removed_df



def convert_target_2_subtype(pos, target_a2r, bridge_r2a, bridge_a2r, subtype_r2a, prot):
    target_raw_pos = convert_pos(pos, target_a2r)
    target_bridge_align_pos = convert_pos(target_raw_pos, bridge_r2a)
    subtype_bridge_raw_pos = convert_pos(target_bridge_align_pos, bridge_a2r)
    subtype_align_pos = convert_pos(subtype_bridge_raw_pos, subtype_r2a)

    return subtype_align_pos


def convert_pos(pos, conv_map):
    # check if position is in map
    if pos in conv_map.keys():
        return conv_map[pos]
    # otherwise, return nearest value
    else:
        return conv_map[min(list(conv_map.keys()), key=lambda x:abs(x-pos))]
    
        


def extract_align_map(prot_fasta_file, rev=False, read_fasta=False):
    prot_align_dict = dict()
    prot_fasta_dict = dict(pd.read_csv(prot_fasta_file, sep='\t').values)
    # loop through each protein
    for prot, fasta in prot_fasta_dict.items():
        # extract coordinate map
        if not read_fasta:
            prot_align_dict[prot] = ft.alignCoordMap(fasta, rev=rev)
        # otherwise, extract fasta dict
        else:
            prot_align_dict[prot] = ft.read_fasta_dict(fasta)

    return prot_align_dict



def create_filepath_map(batch_map_filepath):
    batch_map = pd.read_csv(batch_map_filepath, sep='\t')

    # make sure all entries exist
    for i, row in batch_map.iterrows():
        if not os.path.exists(row["TargetMapFile"]):
            raise FileNotFoundError(f"{row['TargetMapFile']} does not exist.")
        elif not os.path.exists(row["BridgeMapFile"]):
            raise FileNotFoundError(f"{row['BridgeMapFile']} does not exist.")
        elif not os.path.exists(row["SubtypeDir"]):
            raise FileNotFoundError(f"{row['SubtypeDir']} does not exist.")
        elif not os.path.exists(row["EpitopeFile"]):
            raise FileNotFoundError(f"{row['EpitopeFile']} does not exist.")

    return batch_map


def make_dir(path, new):
    dir_name = os.path.join(path, new)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    return dir_name


if __name__ == "__main__":
    main()