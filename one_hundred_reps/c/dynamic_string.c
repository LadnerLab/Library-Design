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
void ds_check_for_resize( dynamic_string_t* input, unsigned int input_len )
{
    int new_capacity;
    int add_length = input_len;
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
    return strlen( input );
}

void ds_init( dynamic_string_t* input )
{
    input->capacity = DEFAULT_LENGTH;
    input->data = calloc( DEFAULT_LENGTH + 1, sizeof( char ) );
    input->size = 0;
}


void ds_add( dynamic_string_t* input, char string[] )
{
    int size = input->size;
    int input_length = strlen( string );
    int new_size = size + input_length;

    ds_check_for_resize( input, input_length );
    strcat( input->data, string );
    input->size = new_size;
}

void ds_clear( dynamic_string_t* input )
{
    free( input->data );
    free( input );
}



