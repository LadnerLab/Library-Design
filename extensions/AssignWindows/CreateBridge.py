#!/usr/bin/env python3
import glob
import os
import fastatools as ft
import argparse
import pandas as pd
import subprocess

def  main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i', '--bridge_map', required=True)
    parser.add_argument('-o', '--output-dir',  help='Output directory for where to output protein combination files.', required=True)
    
    args = parser.parse_args()

    assert not os.path.exists(args.output_dir), \
        f"Warning: The directory \"{args.output_dir}\" already exists. Please move or delete."

    os.mkdir(args.output_dir)

    # create map from old species dir to new speices dirdirectory (checks for valid filepaths)
    bridge_map = create_bridge_map(args.bridge_map)

    create_bridge_alignments_multiple_species(
        bridge_map=bridge_map,
        output_dir=args.output_dir
        )


def create_bridge_alignments_multiple_species(bridge_map, output_dir):
    # create bridge alignments for each batch
    for old_species_file, new_species_file in bridge_map.items():
        bridge_out_dir = make_dir(output_dir, f"Bridge_{old_species_file.split(os.sep)[-2]}_{new_species_file.split(os.sep)[-2]}")

        create_bridge_alignments_single_species(old_species_file, new_species_file, bridge_out_dir)


def create_bridge_alignments_single_species(old_species_file, new_species_file, output_dir):
    # dataframe for mapping protein to file
    target_files_data = list()

    # returns protein to aligned file path dict
    prot_alignment_dict = dict()

    # create alignment dir
    align_dir = make_dir(output_dir, "aligned")

    # generate dictionaries mapping protein to fasta dict
    old_prot_seqs = extract_fasta_dict( old_species_file )
    new_prot_seqs = extract_fasta_dict( new_species_file )

    # make sure proteins are the same
    assert old_prot_seqs.keys() == new_prot_seqs.keys(), \
        "Protein names do not match for all alignments."

    # loop through each protein
    for prot in old_prot_seqs.keys():
        fasta_output_dir = os.path.join(output_dir, f"{prot}.fasta")
        aligned_fasta_output_dir = os.path.join(align_dir, f"{prot}_aligned.fasta")
        bridge_seqs = dict()

        old_seqs = old_prot_seqs[prot]
        new_seqs = new_prot_seqs[prot]

        # get largest seq in each
        old_largest_seq_name = find_largest_seq(old_seqs)
        new_largest_seq_name = find_largest_seq(new_seqs)

        # create fasta dicts
        bridge_seqs[old_largest_seq_name] = old_seqs[old_largest_seq_name].translate({ord('-'): None})
        bridge_seqs[new_largest_seq_name] = new_seqs[new_largest_seq_name].translate({ord('-'): None})

        # write bridge fasta
        ft.write_fasta_dict(bridge_seqs, fasta_output_dir)

        # create aligned fasta
        align_file(
            in_fasta = fasta_output_dir,
            out_fasta = aligned_fasta_output_dir,
            quiet = True
            )

        # append protein and fasta
        target_files_data.append((prot, aligned_fasta_output_dir))
        
    # save target files fasta
    target_file = pd.DataFrame(target_files_data, columns=["Protein", "Fasta"])
    target_file.to_csv(os.path.join(output_dir, f"target_files.tsv"), sep='\t', index=False)


def create_bridge_map(bridge_map_filepath):
    bridge_map = dict(pd.read_csv(bridge_map_filepath, sep='\t').values)

    # make sure all entries exist
    for old_species_dir, new_species_dir in bridge_map.items():
        if not os.path.exists(old_species_dir):
            raise FileNotFoundError(f"{old_species_dir} does not exist.")
        elif not os.path.exists(new_species_dir):
            raise FileNotFoundError(f"{new_species_dir} does not exist.")

    return bridge_map


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


def extract_fasta_dict(prot_fasta_file, rev=False):
    prot_fasta_dict = dict()
    prot_fasta_file_dict = dict(pd.read_csv(prot_fasta_file, sep='\t').values)
    for prot, fasta_file in prot_fasta_file_dict.items():
        prot_fasta_dict[prot] = ft.read_fasta_dict(fasta_file)

    return prot_fasta_dict


def find_largest_seq(fasta_dict):
    # find the largest in each
    largest_seq_len = 0
    for name, seq in fasta_dict.items():
        unaligned_seq = seq.translate({ord('-'): None})
        if len(unaligned_seq) >= largest_seq_len:
            largest_seq_len = len(unaligned_seq)
            largest_seq_name = name

    return largest_seq_name

if __name__ == "__main__":
    main()