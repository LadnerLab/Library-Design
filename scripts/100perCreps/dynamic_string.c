#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

#include "dynamic_string.h"

#define DEFAULT_LENGTH 256
#define SPACE ' '

/**
 * Checks whether or not a dynamic_string_t object needs to be resized
 * if it does, automatically reallocs the object's string buffer
 * to have DEFAULt_LENGTH more characters
 * 
 * @param input dynamic_string_t input to be tested 
 **/
void ds_check_for_resize( dynamic_string_t* input, char string_to_add[] )
{
    int new_capacity;
    int add_length = strlen( string_to_add );
    char* new_data;

    if( input->capacity <= input->size + add_length + 10 )
        {
            new_capacity = ( input->capacity ) + add_length + DEFAULT_LENGTH;
            new_data = realloc( input->data, new_capacity );

            if( !new_data )
                {
                    printf( "Failure to add to dynamic_string\n" );
                }
            else
                {
                    input->data = new_data;
                }

            input->capacity = new_capacity;
        }
}

int string_length( char* input )
{
    int length = 0;

    while( *( input + length ) )
        {
            length++;
        }
    return length;
}

void ds_init( dynamic_string_t* input )
{
    input->capacity = DEFAULT_LENGTH;
    input->data = calloc( sizeof( char ), DEFAULT_LENGTH );
    input->size = 0;
}


void ds_add( dynamic_string_t* input, char string[] )
{
    int size = input->size;
    int input_length = strlen( string );
    int new_size = size + input_length;


    int index = 0;

    ds_check_for_resize( input, string );
 
    for( index = 0; index < input_length; index++ )
        {
            if( string[ index ] >= SPACE )
                {
                    *( input->data + size + index ) = string[ index ];
                }
        }

    input->data[ new_size ] = '\0';
    input->size = new_size;
}

void ds_clear( dynamic_string_t* input )
{
    free( input->data );
    free( input );
}



