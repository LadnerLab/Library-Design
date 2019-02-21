//
//  table.h
//  codon_sampling
//
//  Created by Paul Altin on 11.05.18.
//  Copyright © 2018 Paul Altin. All rights reserved.
//

#ifndef table_h
#define table_h
#include <stdint.h>


const uint8_t MAX_CODONS = 6;
const uint8_t CODON_SIZE = 3;
const uint8_t A_INDEX            = 0;
const uint8_t C_INDEX            = 1;
const uint8_t G_INDEX            = 2;
const uint8_t T_INDEX             = 3;

static_assert(sizeof(char) == 1, "huh?");


struct codon {
    char c[3];
    double w        = 0.0;
    uint8_t index   = 0;
    uint16_t nucleotides[4] = {0};
};


struct table
{
private:
    codon* map['Z'][MAX_CODONS] = {};
    
public:
    
    table(const char* file) {
        
        // read lines from file
        {
            int n, ln = 0; char aa; codon c;
            FILE* f = fopen(file, "r");
            while ((n = fscanf(f, "%c,%3c,%lf,%hhu\n", &aa, c.c, &c.w, &c.index)) != EOF)
                {
                    ++ln;
                    if (n == 4)
                        {
                            printf("Read %c%c%c with weight %f for %c (index %d)\n", c.c[0], c.c[1], c.c[2], c.w, aa, c.index);
                        }
                    else
                        {
                            printf("Unable to read data from line %d, skipping...\n", ln);
                        }
                
                    codon** ptr = &map[ (int) aa][0];
                    while (*ptr) ++ptr;
                    *ptr = new codon(c);
                }
            fclose(f);
        }
        
        uint8_t aa = 0;
        uint8_t i  = 0;
        uint8_t j  = 0;
        // normalize weights and calculate nucleotide count
        for (aa = (uint8_t) 'A'; aa < (uint8_t) 'Z'; ++aa)
            {
                // calculate total of weights
                double total = 0;
                codon** ptr = &map[aa][0];
                while (*ptr && ptr < &map[aa+1][0])
                    {
                        total += (*ptr++)->w;
                    }
            
                // divide each by total
                ptr = &map[aa][0];
                while (*ptr && ptr < &map[aa+1][0])
                    {
                        (*ptr++)->w /= total;
                    }
            
                // add up nucleotide occurrences
                for (i = 0; i < MAX_CODONS; ++i)
                    {
                        for (j = 0; j < 3; ++j)
                            {
                                if (codon* cod = map[aa][i])
                                    {
                                        switch (cod->c[j])
                                            {
                                            case 'A':  ++cod->nucleotides[ A_INDEX ];  break;
                                            case 'C':  ++cod->nucleotides[ C_INDEX ];  break;
                                            case 'G':  ++cod->nucleotides[ G_INDEX ];  break;
                                            case 'T':  ++cod->nucleotides[ T_INDEX ];  break;
                                            }
                                    }
                            }
                    }
            }
    }

    table( const table &tab )
    {
        uint8_t index       = 0;
        uint8_t inner_index = 0;

        for( index = 0; index < 'Z'; index++ )
            {
                for( inner_index = 0; inner_index < MAX_CODONS; inner_index++ )
                    {
                        if( tab.map[ index ][ inner_index ] )
                            {
                                map[ index ][ inner_index ] = new codon( *(tab.map[ index ][ inner_index ]) );
                            }
                    }
            }
    }

    
    ~table()
    {
        for (int i = 0; i < 'Z'; ++i)
            {
                for (int j = 0; j < MAX_CODONS; ++j)
                    {
                        delete map[i][j];
                    }
            }
    }
    
    codon** operator[](char index)
    {
        return &map[ (unsigned int) index ][ 0 ];
    }
};



#endif /* table_h */
