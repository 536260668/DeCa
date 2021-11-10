//
// Created by 梦想家xixi on 2021/11/10.
//

#ifndef MUTECT2CPP_MASTER_ALIGNMENTUTILS_H
#define MUTECT2CPP_MASTER_ALIGNMENTUTILS_H


#include "cigar/Cigar.h"

class AlignmentUtils {
public:
    static Cigar* consolidateCigar(Cigar* c);
    static bool needsConsolidation(Cigar* c);
};


#endif //MUTECT2CPP_MASTER_ALIGNMENTUTILS_H
