//
// Created by 梦想家xixi on 2021/11/25.
//

#ifndef MUTECT2CPP_MASTER_CIGARUTILS_H
#define MUTECT2CPP_MASTER_CIGARUTILS_H

#include "cigar/Cigar.h"
#include "SWParameters.h"
#include "SmithWatermanAlignment.h"

class CigarUtils {
public:
    static Cigar* calculateCigar(uint8_t* refSeq, int refLength, uint8_t* altSeq, int altLength);
    static const SWParameters NEW_SW_PARAMETERS;
    static Cigar* leftAlignCigarSequentially(Cigar* cigar, uint8_t* refSeq, int refLength, uint8_t* readSeq, int readLength, int refIndex, int readIndex);

private:
    static const int SW_PAD = 10;
    static bool isSWFailure(SmithWatermanAlignment* alignment);
};


#endif //MUTECT2CPP_MASTER_CIGARUTILS_H
