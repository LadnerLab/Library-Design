#!/usr/bin/env python3
import tempfile
import pandas as pd
import argparse
import os
from collections import defaultdict

from findEpitopes import read_check_align_file, process_files_probes, iterative_peptide_finder, generate_out_data, create_line_charts
from LinkAlignments import map_bridge_and_link_single_row

WINDOW_SIZE = 30

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i', '--input-table',  help='Tab delimited file with name, target clusters directory file path, and subtype directory file path', required=True)
    parser.add_argument('--cluster-protein-map', help='Tab delimited file with species name, each clusterID, and its corresponding protein.', required=False)
    parser.add_argument('--window-size', type=int, default=WINDOW_SIZE,  help='Size of AA window to use for identifying core epitopes.', required=False)
    parser.add_argument('--max-zeros', type=int, default=5, help='Maximum number of zero counts a window can contain.', required=False)
    parser.add_argument('--max-overlap', type=int, default=8, help='Maximum AA overlap a window can have with a previously selected window.', required=False)
    parser.add_argument('--peptide-overlap', type=int, default=9, help='Peptide sequence should overlap at least this amount to be included in output data.', required=False)
    parser.add_argument('--peak-overlap-window-size', type=int, default=10, help='Window size around found peaks in which containing peptides will be removed for the next iteration.', required=False)
    parser.add_argument('--include-iter-vis', action="store_true", help='Output each chart given the iteration.', required=False)
    parser.add_argument('-o', '--output-dir', default="find_epitopes_and_link_out", help='Name of directory to output line plots, output protein combination files, and bridge alignments.')
    
    parser.add_argument('-k', '--kmer-size', type=int, default=6, help='Used for assigning sequences to proteins. Kmer size for kmer sets that are going to be intersected.', required=False)
    parser.add_argument('--kmer-ovlp-thresh', type=float, default=0.3, help='Used for assigning sequences to proteins. Minimum kmer overlap that the largest overlapping sequence must have to be assigned to a protein.', required=False)
    parser.add_argument('--max-window-gaps', type=int, default=WINDOW_SIZE/2, required=False)

    args = parser.parse_args()


    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    else:
        print(f"Warning: The directory \"{args.output_dir}\" already exists. Files may be overwritten.\n")

    # extract data from input table
    input_table = create_filepath_map(args.input_table)

    if args.cluster_protein_map:
        # create cluster protein map dict
        cluster_protein_map = defaultdict(dict)
        # read cluster protein map
        cluster_protein_map_df = pd.read_csv(args.cluster_protein_map, sep='\t')
        for i, row in cluster_protein_map_df.iterrows():
            cluster_protein_map[row["Name"]][row["ClusterID"]] = row["Protein"]

    # loop through each row in input table
    for i, row in input_table.iterrows():
        spec_output_dir = make_dir(args.output_dir, row["Name"])

        directory_path = row["TargetClusters"]
        alignment_to_use_dict = read_check_align_file(directory_path)

        # get original align counts and peptide positions 
        alignCountsD, file_2_pep_pos_dict = process_files_probes(
                                                probes_dict=alignment_to_use_dict, 
                                                directory_path=directory_path
                                                )

        windows, last_iter_counts = iterative_peptide_finder(alignment_to_use_dict, directory_path, args.window_size, args.max_zeros, args.max_overlap, 
                                                                    args.peak_overlap_window_size, spec_output_dir, args.include_iter_vis)

        generate_out_data(spec_output_dir, directory_path, last_iter_counts, alignCountsD, file_2_pep_pos_dict, windows, args.window_size, args.peak_overlap_window_size)

        create_line_charts(alignCountsD, windows, spec_output_dir)


        # create temporary target alignments file and use the same for bridge (use clusters = [f"{file.split('_')[-2]} for file in alignCountsD.keys()])
        with tempfile.TemporaryDirectory() as tmpdirname:

            target_files_tsv = os.path.join(tmpdirname, "target_files.tsv")

            with open(target_files_tsv, 'w') as target_files:

                target_files.write("Protein\tFasta")

                # generate target files tsv
                for file in alignCountsD.keys():
                        target_files.write(f"\n{file.split('_')[-2]}\t{os.path.join(directory_path, file)}")


            epitope_positions_file = os.path.join(spec_output_dir, "peptide_seq_data.tsv")

            map_bridge_and_link_single_row(
                kmer_size=args.kmer_size,
                kmer_ovlp_thresh=args.kmer_ovlp_thresh,
                max_window_gaps=args.max_window_gaps,
                output_dir=args.output_dir,
                output_name=row["Name"],
                target_alignments_file=target_files_tsv,
                subtype_dir=row["SubtypeDir"],
                bridge_alignments_file=target_files_tsv,
                epitope_positions_file=epitope_positions_file
                )

        # overwrite output tables cluster id
        if args.cluster_protein_map:
            add_protein_name(epitope_positions_file, cluster_protein_map, row["Name"])
            add_protein_name(os.path.join(spec_output_dir, "new_peptide_seq_data.tsv"), cluster_protein_map, row["Name"])
            add_protein_name(os.path.join(spec_output_dir, "removed_peptides.tsv"), cluster_protein_map, row["Name"])

        
def add_protein_name(epitope_positions_file, cluster_protein_map, spec_name):
    temp_rename_epitope_pos_df = pd.read_csv(epitope_positions_file, sep='\t')
    # add row
    temp_rename_epitope_pos_df.insert(loc=0, column="Protein", value=[None]*len(temp_rename_epitope_pos_df.index))
    for i, pep_row in temp_rename_epitope_pos_df.iterrows():
        if pep_row["ClusterID"] in cluster_protein_map[spec_name].keys():
            temp_rename_epitope_pos_df.at[i, "Protein"] = cluster_protein_map[spec_name][pep_row["ClusterID"]]

    temp_rename_epitope_pos_df.to_csv(epitope_positions_file, sep='\t', index=False)


def create_filepath_map(batch_map_filepath):
    batch_map = pd.read_csv(batch_map_filepath, sep='\t')

    # make sure all entries exist
    for i, row in batch_map.iterrows():
        if not os.path.exists(row["TargetClusters"]):
            raise FileNotFoundError(f"{row['TargetClusters']} does not exist.")
        elif not os.path.exists(row["SubtypeDir"]):
            raise FileNotFoundError(f"{row['SubtypeDir']} does not exist.")

    return batch_map


def make_dir(path, new):
    dir_name = os.path.join(path, new)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    return dir_name

if __name__ == "__main__":
    main()