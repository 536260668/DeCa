//
// Created by 梦想家xixi on 2021/11/25.
//

#ifndef MUTECT2CPP_MASTER_CIGARUTILS_H
#define MUTECT2CPP_MASTER_CIGARUTILS_H

#include "cigar/Cigar.h"

class CigarUtils {
public:
    static Cigar* calculateCigar(uint8_t* refSeq, int refLength, uint8_t* altSeq, int altLength);
};


#endif //MUTECT2CPP_MASTER_CIGARUTILS_H
