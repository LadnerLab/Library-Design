//
//  xoroshiro.h
//  codon_sampling
//
//  Created by Paul Altin on 10.05.18.
//  Copyright © 2018 Paul Altin. All rights reserved.
//

#ifndef xoroshiro_h
#define xoroshiro_h

#include <random>


// this is xoroshiro128+ 1.0, by David Blackman and Sebastiano Vigna
// stolen from http://xoroshiro.di.unimi.it/xoroshiro128plus.c

namespace xoroshiro
{
    static uint64_t s[2] = { 0x41, 0x29837592 };
    
    inline void seedrandom() {
        std::random_device rd;
        uint64_t seed = rd();
        seed = (seed << 32) | rd();
        s[0] = (1181783497276652981ULL * seed);
        s[1] = (1181783497276652981ULL * s[0]);
    }
    
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
    
    inline uint64_t next(void) {
        const uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        const uint64_t result = s0 + s1;
        s1 ^= s0;
        s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
        s[1] = rotl(s1, 36); // c
        return result;
    }
    
    inline double uniform() {
        return (1.0/18446744073709551616.0) * next();
    }
}
    

#endif /* xoroshiro_h */
