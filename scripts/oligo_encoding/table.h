//
//  table.h
//  codon_sampling
//
//  Created by Paul Altin on 11.05.18.
//  Copyright © 2018 Paul Altin. All rights reserved.
//

#ifndef table_h
#define table_h


const size_t MAX_CODONS = 6;

static_assert(sizeof(char) == 1, "huh?");


struct codon {
    char c[3];
    double w = 0.0;
    int index = 0;
    int nucleotides[4] = {0};
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
            while ((n = fscanf(f, "%c,%3c,%lf,%d\n", &aa, c.c, &c.w, &c.index)) != EOF)
            {
                ++ln;
                if (n == 4)
                    printf("Read %c%c%c with weight %f for %c (index %d)\n", c.c[0], c.c[1], c.c[2], c.w, aa, c.index);
                else
                    printf("Unable to read data from line %d, skipping...\n", ln);
                
                codon** ptr = &map[aa][0];
                while (*ptr) ++ptr;
                *ptr = new codon(c);
            }
            fclose(f);
        }
        
        // normalize weights and calculate nucleotide count
        for (char aa = 'A'; aa < 'Z'; ++aa)
        {
            // calculate total of weights
            double total = 0;
            codon** ptr = &map[aa][0];
            while (*ptr && ptr < &map[aa+1][0])
                total += (*ptr++)->w;
            
            // divide each by total
            ptr = &map[aa][0];
            while (*ptr && ptr < &map[aa+1][0])
                (*ptr++)->w /= total;
            
            // add up nucleotide occurrences
            for (int i = 0; i < MAX_CODONS; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (codon* cod = map[aa][i]) {
                        switch (cod->c[j]) {
                            case 'A':  ++cod->nucleotides[0];  break;
                            case 'C':  ++cod->nucleotides[1];  break;
                            case 'G':  ++cod->nucleotides[2];  break;
                            case 'T':  ++cod->nucleotides[3];  break;
                        }
                    }
                }
            }
        }
    }
    
    ~table() {
        for (int i = 0; i < 'Z'; ++i)
            for (int j = 0; j < MAX_CODONS; ++j)
                delete map[i][j];
    }
    
    codon** operator[](char index) {
        return &map[index][0];
    }
};



#endif /* table_h */
