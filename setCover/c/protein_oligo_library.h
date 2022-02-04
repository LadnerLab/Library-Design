#ifndef PROTEIN_OLIGO_H_INCLUDED
#define PROTEIN_OLIGO_H_INCLUDED

#include "dynamic_string.h"
#include "array_list.h"
#include "hash_table.h"
#include "set.h"

typedef struct sequence_t
{
    char* name;
    dynamic_string_t* sequence;
} sequence_t;

typedef struct subset_data_t
{
    unsigned int start;
    unsigned int end;
} subset_data_t;

typedef struct blosum_data_t
{
    char* letter_data;
    hash_table_t* blosum_table;
    
} blosum_data_t;



/**
 * Reads all of the fastas in given file
 * @param file_to_read string name of file to read
 * @return array of pointers to sequences containing sequence data
 *         found in file
 **/
void read_sequences( FILE* file_to_read, sequence_t** in_sequence );

/**
 * Counts the amount of fastas in a fasta file provided
 * 
 * Note: Resets file pointer
 * @param file_to_read fasta file to be opened and read
 * @returns int number of sequences in file, or -1 if
 *          file could not be opened
**/
int count_seqs_in_file( FILE* data_file );

/**
 * Reads a line from a file into a new dynamic_string_t object
 * @param stream open file pointer to read from
 * @returns dynamic_string_t representation of the line
 **/
int get_a_line( FILE* stream, dynamic_string_t* to_read );

/**
 * Counts the number of specified character in a string
 * @param string_in pointer to first character of a string
 * @param to_find character to count in the string
 * @return integer number of to_find character found in string
 **/
int count_char_in_string( char* string_in, char to_find );


/**
 * Tests for character in a string
 * @param string_in pointer to first character in a string
 * @param to_find character value to test for
 * @returns boolean result of test
 **/
int char_in_string( char* string_in, char to_find );

/**
 * Calculates the percentage of specified charater makes up string
 * @param string_in string to be tested
 * @param test_char character to test for in string
 * @return floating point percentage of character in string
 **/
float percent_char_in_string( char* string_in, char test_char );

/**
 * Writes fastas to an output file, from an array of sequence_ts.
 * @param in_seqs array of pointers to sequences to write
 * @param num_seqs the number of sequences to write to the file
 * @param output string filename to write to
 **/
void write_fastas( sequence_t** in_seqs, int num_seqs, char* output );

/**
 * Removes a character from a string in place
 * Note: String passed into function must be declared as a character array,
 * not a character pointer.
 * @param test_string String from which to remove the character
 * @param to_remove character to remove from string 
 **/
void remove_char_from_string( char* test_string, char to_remove );

/**
 * Determines whether a given sequence is or is not valid
 * @param sequence String sequence to be tested
 * @param minimum length of non dash characters needed to be present
          in order for the string to count
 * @param percent_valid Percentage of characters needed to be non-dash in order
          for the string to count
**/
int is_valid_sequence( char* sequence, int min_length, float percent_valid );
/**
 * Appends all valid xmers within a sequence to a hashtable
 * @param in_hash pointer to hash_table to add the valid xmers to
 * @param in_seq pointer to string to create a subset of 
 * @param window_size integer number of characters to capture with each iteration
 * @param step_size integer number of characters to move over after each iteration
 * @returns pointer to hash table containing all of the subsets of the sequence, 
 *          as key, and an array list of subset_data_t as key containing start/end
 **/ 
hash_table_t* create_xmers_with_locs( hash_table_t* in_hash, char* name,
                                      char* in_seq,
                                      int window_size, int step_size );

/**
 * Break a ymer down into into the unique locations of its xmers
 * @param in_ymer_name string name of the ymer input to method
 * @param in_ymer string ymer
 * @param in_xmer_table pointer to xmer table containing xmers to search
 * @param window_size integer size of each xmer
 * @param step_size integer amount to move over after each ymer capture
 * @param permute boolean option to permute xmer functional groups
 * @param blosum_data pointer to blosum_data_t containing the data found from a blosum
          matrix file
 *         Note: if the -b flag is not used, this will be NULL
 * @param blosum_cutoff integer cutoff score that means a change
 *        will be made

 * @returns set of strings containing the locations of in_ymer's xmers
 **/
set_t* component_xmer_locs( char* in_ymer_name, char* in_ymer,
                            set_t* out_ymer,
                            hash_table_t* in_xmer_table,
                            int window_size, int step_size,
                            blosum_data_t* blosum_data,
                            int blosum_cutoff,
                            int permute
                          );
/**
 * Break a ymer down into its xmers, the value of each xmer is NULL
 * @param in_hash pointer to hash_table to add the valid xmers to
 * @param in_seq pointer to string to create a subset of 
 * @param window_size integer number of characters to capture with each iteration
 * @param step_size integer number of characters to move over after each iteration
 * @param permute boolean option whether or not xmer permutation should be done
 * @param blosum_data pointer to blosum_data_t containing the data found from a blosum
 *         matrix file
 *         Note: if the -b flag is not used, this will be NULL
 * @param blosum_cutoff integer cutoff score that means a change
 *        will be made

 * @returns pointer to hash table containing all of the subsets of the sequence, 
 *          as key, and an array list of subset_data_t as key containing start/end
 **/
hash_table_t* subset_lists( hash_table_t* in_hash,
                            char* in_seq,
                            int window_size, int step_size,
                            blosum_data_t* blosum_data,
                            int blosum_cutoff,
                            int permute
                          );


/**
 *  Frees all of the pointers found in in_data array_list
 *  @param in_data array_list of pointers to be free'd
 **/
void free_data( array_list_t* in_data );

/**
 * allocates num_bytes of memory, stores the pointer in an array list
 * and returns the pointer given by malloc
 * @param data pointer to array_list to hold pointers
 * @param num_bytes number of bytes to allocate
 * @returns integer boolean success of operation 
**/
void *malloc_track( array_list_t* data, int num_bytes );

/**
 * Calculates the number of subsequences created, using
 * length of the sequence and the stepsize
 * @param length integer length of the sequence to calculate
 * 		  number of subsequences for
 * @param window_size integer number of sequences to capture per iteration
 * @returns integer number of subsequences created
 * */
int calc_num_subseqs( int length, int window_size );


/**
 * Creates each 1-permutation of an kmer, based on the functional groupings of 
 * amino acids.
 * @param str_to_change kmer to change
 * @param permutations pointer to array_list_t in which to store a kmer's permutations
 * @param blosum_data pointer to blosum_data_t containing the data found from a blosum
 *         matrix file
 *         Note: if the -b flag is not used, this will be NULL
 * @param blosum_cutoff integer cutoff score that means a change
 *        will be made
 **/
void permute_xmer_functional_groups( char* str_to_change,
                                     array_list_t* permutations,
                                     blosum_data_t* blosum_data,
                                     int blosum_cutoff
                                   );

/**
 * Parses a blosum file containing a blosum substitution
 * matrix
 * @param file_name string name of file to open and parse
 * @returns blosum_data_t* pointer to struct containing
 *          an array of character amino acids found in the file,
 *          and integer array of distances between these amino acids
 **/
blosum_data_t* parse_blosum_file( FILE* file_name );

#endif
