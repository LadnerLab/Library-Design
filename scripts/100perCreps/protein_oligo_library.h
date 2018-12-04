#ifndef PROTEIN_OLIGO_H_INCLUDED
#define PROTEIN_OLIGO_H_INCLUDED

#include "dynamic_string.h"
#include "array_list.h"

typedef struct sequence_t
{
    char* name;
    dynamic_string_t* sequence;
} sequence_t;


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
 * Calculates the number of subsequences created, using
 * length of the sequence and the stepsize
 * @param length integer length of the sequence to calculate
 * 		  number of subsequences for
 * @param window_size integer number of sequences to capture per iteration
 * @returns integer number of subsequences created
 * */
int calc_num_subseqs( int length, int window_size );


#endif
