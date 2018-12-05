#ifndef HASHTTABLE_HH_INCLUDED
#define HASHTTABLE_HH_INCLUDED

#include <stdint.h>
#define ITEM_NOT_FOUND -1

typedef struct HT_Entry
{
    char* key;
    void* value;
    struct HT_Entry* next;
    struct HT_Entry* prev;
} HT_Entry;

typedef struct hash_table_t
{
    HT_Entry** table_data; 
    uint32_t size;
    uint32_t capacity;
} hash_table_t;


/**
 * Initializes a HastTable struct. Allocates memory
 * to store size number of HT_Entry structs.
 * @param table pointer to hash_table_t to init
 * @param size integer size of hash_table_t
 **/
void ht_init( hash_table_t* table, int size );

/**
 * Clears a hash_table_t struct. Frees the memory
 * allocated for each HT_entry
 * @param table hash_table_t object to free
 **/
void ht_clear( hash_table_t* table ); 
/**
 * Generates the hash for use by the table.
 * @param input string to calculate hash for 
 * @return integer value representing hash of object
 **/
uint32_t generate_hash( const void *key,  int len, uint32_t seed );

/**
 * Adds an entry to the hashtable
 * Note: Uses linked list to resolve collisions
 * Note: If a duplicate key is added, value at that key is updated 
 *       to provided value
 * @param table pointer to hash_table_t to add data to
 * @param to_add string value to add as key to hash_table
 * @param add_val value to be added as value of the key
 * @returns integer value representing success of addition to table
 **/
int ht_add( hash_table_t* table, char* to_add, void* add_val );


/**
 * Finds the index of an item in the hash table
 * Note: Uses linked list to resolve collisions
 * @param table pointer to hash_table_t whose capacity is used to 
 *        calculate initial index
 * @param in_key string key value to be searched for in the table
 * @returns pointer to HT_Entry found in list, or NULL if item was not found
 **/
HT_Entry* find_item( hash_table_t* table, char* in_key );

/**
 * Finds an item within a hash table using provided key.
 * Note: Uses find_item
 * @param table pointer to hash_table_t to search
 * @param in_key String key value to use to search
 * @returns pointer to the integer found with the key in_key
 *          or NULL if the value was not found
 **/
void *ht_find( hash_table_t* table, char* in_key );


/**
 * Removes and frees data found in a hash table
 * Note: Uses find_item
 * @param table pointer to hash_table_t to delete from
 * @param in_key String key value to search
 * @returns pointer to value of node that was deleted
 **/
void* ht_delete( hash_table_t* table, char* in_key );

/** 
 * Gets all of the HT_Entry items, stores them in an array
 * @param input pointer to hash_table_t to get the items of
 * @returns pointer to array of pointers to HT_Entry items
 **/ 
HT_Entry *ht_get_items( hash_table_t* input );
HT_Entry *ht_get_items_no_malloc( hash_table_t* input, HT_Entry *data );
#endif

