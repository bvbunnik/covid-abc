#define __STDC_LIMIT_MACROS
#include <iostream>

#include "SFMT/SFMT.h"

#include "../include/u01.hpp"

uint64_t u01n;
uint32_t *rnd_array=NULL;
sfmt_t sfmt;

uint32_t *skew_seed(bool random) {
    uint32_t *seed=new uint32_t[SKEW_SEEDS];
    time_t tick;
    uint j=0;
    for (uint i=0;i<SKEW_SEEDS;i++) {
        if (random) {
            tick=clock();
            for(j=0;(uint)clock()==tick;j++)
                    ;
        }

            seed[i]= j%256;// only use lower 8 bits..
    }
    return seed;
}

float u01(char c) {

    float f;

    repeat:
        u01n++;

    if (u01n%RND_BUFFER==1) {
        if (rnd_array==NULL) {

            if (np==0) {
                std::cerr<<"+------------------------------------------------------------------------------+\n";
                std::cerr<<"|             Invoking SIMD-oriented Fast Mersenne Twister (SFMT)              |\n";
                std::cerr<<"| Copyright 2006,2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima University  |\n";
                std::cerr<<"| Copyright 2012 Mutsuo Saito, Makoto Matsumoto, Hiroshima University and      |\n";
                std::cerr<<"|                The University of Tokyo. All rights reserved.                 |\n";
                std::cerr<<"+------------------------------------------------------------------------------+"<<std::endl;
            }

            rnd_array=new uint32_t[RND_BUFFER];
            uint32_t *seeds=skew_seed(RANDOM_SEED);
            sfmt_init_by_array(&sfmt, seeds, SKEW_SEEDS);
            delete [] seeds;
        }
        sfmt_fill_array32(&sfmt, rnd_array, RND_BUFFER);
        u01n=1;
    }

    f=rnd_array[u01n-1]/(float)UINT32_MAX;

    switch (c) {
        case LEFT_OPEN:
            if (f==0)
                goto repeat;
            else
                return f;
        case RIGHT_OPEN:
            if (f==1)
                goto repeat;
            else
                return f;
        case BOTH_OPEN:
            if (f==0 || f==1)
                goto repeat;
            else
                return f;
        default:
            return f;
    }
}

float u01() { //this is a function overload... will give in [0,1] interval.
    return u01(0);
}
