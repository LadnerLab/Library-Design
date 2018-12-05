#ifndef DYNAMIC_STRING_H_INCLUDED
#define DYNAMIC_STRING_H_INCLUDED
#define DEFAULT_LENGTH 256

typedef struct dynamic_string_t
{
    int size;
    int capacity;
    char* data; 
} dynamic_string_t;

/**
 * Calculates the length of a string
 * @param input character pointer to first character of string
 * @return integer length of string
 *
 **/
int string_length( char* input );

/**
 * Checks whether or not a dynamic_string_t object needs to be resized
 * if it does, automatically reallocs the object's string buffer
 * to have DEFAULt_LENGTH more characters
 * 
 * @param input dynamic_string_t input to be tested 
 **/
void ds_check_for_resize( dynamic_string_t* input, char string_to_add[] );


/**
 * Initializes input dynamic_string_t struct instance
 * Note: mallocs memory on the heap for object's data
 * @param input pointer to dynamic_string_t instance to be initialized
 **/
void ds_init( dynamic_string_t* input );

/**
 * Adds new string to a dynamicString object
 * Note: Uses check_for_resize method to see if the object should
 * be given more memory
 * @param input dynamic_string_t object to be added to
 * @param string[] string array to append to input
 **/
void ds_add( dynamic_string_t* input, char string[] );

/**
 * clears a dynamic_string_t struct instance
 * Note: Frees both the object itself and the data
 * it contains
 * @param input pointer to dynamic_string_t instance to be freed
 **/
void ds_clear( dynamic_string_t* input );


#endif
