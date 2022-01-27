#define _GNU_SOURCE
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "protein_oligo_library.h"
#include "dynamic_string.h"
#include "array_list.h"

#define LINE_SIZE 512
#define DASH_CHAR '-'
#define SPACE ' '
#define DATA_NOT_FOUND -99

// ============== Local Function Prototypes ==================== // 

int count_letters( char* in_str );



// ============================================================= // 

int num_digits_in_int( int input )
{
    char int_as_str[ LINE_SIZE ];
    return sprintf( int_as_str, "%d", input );
}

int count_letters( char* in_str )
{
    int char_count = 0;
    int index = 0;

    while( in_str[ index ] )
        {
            if( in_str[ index ] >= 'A' && in_str[ index ] <= 'Z' )
                {
                    char_count++;
                }
            index++;
        }

    return char_count;
}

void read_sequences( FILE* file_to_read, sequence_t** in_sequence )
{
    int has_line;
    int index = 0;

    dynamic_string_t* line = (dynamic_string_t*) malloc( sizeof( dynamic_string_t ) );
    dynamic_string_t* sequence = NULL;

    ds_init( line );
    has_line = get_a_line( file_to_read, line );

    while( has_line )
        {

            // remove newline-character if the line is not empty
            if( strchr( line->data, '\n' ) != NULL
                && line->size )
                {
                    line->data[ --line->size ] = '\0';
                }

            if( line->data[ 0 ] == '>' )
                {
                    *( in_sequence + index ) = malloc( sizeof( sequence_t ) );
                    sequence = malloc( sizeof( dynamic_string_t ) );
                    ds_init( sequence );

                    in_sequence[ index ]->name = line->data;
                    in_sequence[ index ]->sequence = sequence; 
                    in_sequence[ index ]->collapsed = 0;
                    index++;
                }
            else
                {
                    ds_add( sequence, line->data );
                }
            ds_init( line );
            has_line = get_a_line( file_to_read, line );

            // remove newline-character, update size to reflect the removal of
            // char
            if( line->size )
                {
                    line->data[ --line->size ] = '\0';
                }
        }

    ds_clear( line );
}



void write_fastas( sequence_t** in_seqs, int num_seqs, char* output )
{
    FILE* out_file = fopen( output, "w+" );
    int index;

    if( !out_file )
        {
            printf( "Unable to open file %s for output." , output );
        }

    for( index = 0; index < num_seqs; index++ )
        {
            fprintf( out_file, "%s\n", in_seqs[ index ]->name );

            fprintf( out_file, "%s\n", in_seqs[ index ]->sequence->data );
        }

    fclose( out_file );
}


int count_seqs_in_file( FILE* data_file )
{

    int counter = 0;
    int current_char = 0;

    if( !data_file )
        {
            return -1; 
        }

    rewind( data_file );

    current_char = fgetc( data_file );
    while( current_char != EOF )
        {
            if( (char) current_char == '>' )
                {
                    counter++;
                }
            current_char =  fgetc( data_file );
        }
    rewind( data_file );
    return counter;
}


int get_a_line( FILE* stream, dynamic_string_t* to_read )
{
    char *current_char = NULL;
    size_t size = 0;
    bool return_val = false;

    if( getline( &current_char, &size, stream ) != EOF ) 
        {
            ds_add( to_read, current_char );

            return_val = true;
        }

    free( current_char );
    return return_val;

}


int count_char_in_string( char* string_in, char to_find )
{
    int count = 0;
    int index = 0;

    while( *(string_in + index ) )
        {
            if( *( string_in + index ) == to_find )
                {
                    count++;
                }
            index++;
        }

    return count;
}

int char_in_string( char* string_in, char to_find )
{
    int index = 0;
    while( *( string_in + index ) )
        {
            if ( *(string_in + index ) == to_find )
                {
                    return true;
                }
            index++;
        }
    return false;
}




