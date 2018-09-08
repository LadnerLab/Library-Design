#!/usr/bin/env python3

import os
import sys
import matplotlib.pyplot as plt
import argparse 

import protein_oligo_library as oligo

def main():
    parser = argparse.ArgumentParser(
             description = "Plot ratio of xmers covered to percent xmer coverage"
    )

    parser.add_argument( '--kmer_suffix', type = str )
    parser.add_argument( '--alignment_suffix', type = str )
    parser.add_argument( '--spanning_dir', type = str )
    parser.add_argument( '--kmer_dir', type = str )
    parser.add_argument( '--ref_dir', type = str )

    args = parser.parse_args()

    spanning_dir = args.spanning_dir
    kmer_dir = args.kmer_dir
    aligned_suffix = args.alignment_suffix
    ref_dir = args.ref_dir

    alignment_size_dict = {}
    alignment_count_dict = {}

    kmer_size_dict = {}
    kmer_count_dict = {}

    ref_dict = {}


    for current_cluster in os.listdir( spanning_dir ):
        names, sequences = oligo.read_fasta_lists( spanning_dir + '/' + current_cluster )

        alignment_size_dict[ current_cluster.split( aligned_suffix )[0] ] = set()
        alignment_count_dict[ current_cluster.split( aligned_suffix )[0] ] = len( names )

        for current_seq in sequences:
            alignment_size_dict[ current_cluster.split(aligned_suffix )[0] ] |= oligo.subset_lists_iter( current_seq, 10, 1 )


    for current_cluster in os.listdir( kmer_dir ):
        names, sequences = oligo.read_fasta_lists( kmer_dir + '/' + current_cluster )

        kmer_size_dict[ current_cluster.split('_out_R_1' )[0] ] = set()

        kmer_count_dict[ current_cluster.split('_out_R_1' )[0] ] = len( names )

        for current_seq in sequences:
            kmer_size_dict[ current_cluster.split('_out_R_1' )[0] ] |= oligo.subset_lists_iter( current_seq, 10, 1 )

    for current_cluster in os.listdir( ref_dir ):
        names, sequences = oligo.read_fasta_lists( ref_dir + '/' + current_cluster )

        ref_dict[ current_cluster ] = set()

        for current_seq in sequences:
            ref_dict[ current_cluster ] |= oligo.subset_lists_iter( current_seq, 10, 1 )
        

    xaxis_vals = list()
    yaxis_vals = list()

    cluster_names = list( alignment_size_dict.keys() )

    for current_clust in cluster_names:
        if current_clust in kmer_size_dict.keys():
            if alignment_count_dict[ current_clust ] > 0:
                yaxis_vals.append( kmer_count_dict[ current_clust ] / alignment_count_dict[ current_clust ] )

                percent_kmer_cov = len( kmer_size_dict[ current_clust ] ) / len( ref_dict[ current_clust ] )
                percent_alignment_cov = len( alignment_size_dict[ current_clust ] ) / len( ref_dict[ current_clust ] )


                xaxis_vals.append( percent_kmer_cov / percent_alignment_cov )


    ax = plt.subplot()
    ax.scatter( xaxis_vals, yaxis_vals )
    plt.xlabel( "Percent coverage (kmer / alignment )" )
    plt.ylabel( "Oligo counts (kmer / alignment)" )
    # plt.plot( [ item for item in range( 12000 ) ] )
    plt.ylim( ymax = 12 )
    plt.xlim( xmin = 0, xmax = 12 )

    plt.title( "Gap Spanning" )
    plt.show()






if __name__ == '__main__':
    main()
