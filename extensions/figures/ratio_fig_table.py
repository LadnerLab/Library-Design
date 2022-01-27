#!/usr/bin/env python3

import os
import sys
import argparse 

try:
    import matplotlib.pyplot as plt
except:
    pass

def main():
    parser = argparse.ArgumentParser(
             description = "Plot ratio of xmers covered to percent xmer coverage"
    )

    parser.add_argument( '--kmer_suffix', type = str )
    parser.add_argument( '--alignment_suffix', type = str )
    parser.add_argument( '--spanning_dir', type = str )
    parser.add_argument( '--kmer_dir', type = str )
    parser.add_argument( '--ref_dir', type = str )
    parser.add_argument( '--out', type = str )
    parser.add_argument( '--plot', action = "store_true", default = False )

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
    ref_count_dict = {}
    ref_length_dict = {}
    ref_19mer_len = {}


    for current_cluster in os.listdir( spanning_dir ):
        names, sequences = read_fasta_lists( spanning_dir + '/' + current_cluster )

        alignment_size_dict[ current_cluster.split( aligned_suffix )[0] ] = set()
        alignment_count_dict[ current_cluster.split( aligned_suffix )[0] ] = len( names )

        for current_seq in sequences:
            alignment_size_dict[ current_cluster.split( aligned_suffix )[0] ] |= subset_lists_iter( current_seq, 10, 1 )


    for current_cluster in os.listdir( kmer_dir ):
        names, sequences = read_fasta_lists( kmer_dir + '/' + current_cluster )

        kmer_size_dict[ current_cluster.split('_out_R_1' )[0] ] = set()

        kmer_count_dict[ current_cluster.split('_out_R_1' )[0] ] = len( names )


        for current_seq in sequences:
            kmer_size_dict[ current_cluster.split('_out_R_1' )[0] ] |= subset_lists_iter( current_seq, 10, 1 )

    for current_cluster in os.listdir( ref_dir ):
        names, sequences = read_fasta_lists( ref_dir + '/' + current_cluster )

        ref_count_dict[current_cluster] = len(names)
        ref_dict[ current_cluster ] = set()
        ref_length_dict[ current_cluster ] = 0
        ref_19mer_len[ current_cluster ] = set()

        for current_seq in sequences:
            ref_dict[ current_cluster ] |= subset_lists_iter( current_seq, 10, 1 )
            ref_length_dict[ current_cluster ] += len( current_seq )
            ref_19mer_len[ current_cluster ] |= subset_lists_iter( current_seq, 19, 1  )

    xaxis_vals = list()
    yaxis_vals = list()

    cluster_names = list( alignment_size_dict.keys() )
    kmer_cluster_names = list( kmer_size_dict.keys() )

    #Report clusters not represented in both libraries being compared
    align_only = set(cluster_names).difference(set(kmer_cluster_names))
    if align_only: print ("Clusters absent from kmer library: %s" % (",".join(align_only)))
    kmer_only = set(kmer_cluster_names).difference(set(cluster_names))  
    if kmer_only: print ("Clusters absent from aligned library: %s" % (",".join(kmer_only)))
    
    fout = open(args.out, "w")
    fout.write("Cluster\tNumSeqs\tKmerOligos\tAlignOligos\tOligosRatios\tRefEpitopes\tKmerEpitopes\tAlignEpitopes\tKmerPropEpi\tAlignPropEpif\tPropEpiRatio\tNumUniqYmersInClust\tTotalLenSeqsInCluster\n")
    for current_clust in cluster_names:
        if current_clust in kmer_size_dict.keys():
            if alignment_count_dict[ current_clust ] > 0:

                current_clust_num_kmers = kmer_count_dict[ current_clust ]
                current_clust_num_alignment_kmers = alignment_count_dict[ current_clust ]

                yaxis_vals.append( kmer_count_dict[ current_clust ] / alignment_count_dict[ current_clust ] )

                percent_kmer_cov = len( kmer_size_dict[ current_clust ] ) / len( ref_dict[ current_clust ] )
                percent_alignment_cov = len( alignment_size_dict[ current_clust ] ) / len( ref_dict[ current_clust ] )

                xaxis_vals.append( percent_kmer_cov / percent_alignment_cov )
                
                #Write values to output file
                fout.write("%s\t%d\t%d\t%d\t%.5f\t%d\t%d\t%d\t%.5f\t%.5f\t%.5f\t%d\t%d\n" %
                              (
                                  current_clust, ref_count_dict[current_clust],
                                  kmer_count_dict[current_clust], alignment_count_dict[current_clust], 
                                  kmer_count_dict[ current_clust ] / alignment_count_dict[ current_clust ],
                                  len(ref_dict[current_clust]), len(kmer_size_dict[current_clust]),
                                  len(alignment_size_dict[current_clust]),
                                  percent_kmer_cov, percent_alignment_cov,
                                  percent_kmer_cov/percent_alignment_cov,
                                  len( ref_19mer_len[ current_clust ] ),
                                  ref_length_dict[ current_clust ]
                              )
                )
    
    fout.close()
    
    if args.plot:
        ax = plt.subplot()
        ax.scatter( xaxis_vals, yaxis_vals )
        plt.xlabel( "Percent coverage (kmer / alignment )" )
        plt.ylabel( "Oligo counts (kmer / alignment)" )
        # plt.plot( [ item for item in range( 12000 ) ] )
        plt.ylim( ymax = 12 )
        plt.xlim( xmin = 0, xmax = 12 )
        
        plt.title( "Gap Spanning" )
        plt.show()



def subset_lists_iter(sequence, window_size, step_size ):
    xmer_set = set()

    start = 0
    end = window_size

    while end <= len( sequence ):
        xmer = sequence[start:end]
        if 'X' not in xmer:
            xmer_set.add(xmer)
        start += step_size
        end = start + step_size + window_size - 1

    return xmer_set

def read_fasta_lists( file_to_read ):
    """
       Reads a list of fastas from file_to_read
    
       Returns:
        names- a list of names of the sequences found in the fasta file
        sequences- a list of the sequences found in the fasta file
    """

    file_in = open( file_to_read, 'r' )
    count = 0

    names = []
    sequences = []
    current_sequence = ''

    for line in file_in:
        line = line.strip()
        if line and line[ 0 ] == '>':
            count += 1
            names.append( line[ 1: ] )
            if count > 1:
                sequences.append( current_sequence )
            current_sequence = ''

        else:
            current_sequence += line

    sequences.append( current_sequence )
    file_in.close()

    return names, sequences


if __name__ == '__main__':
    main()
