//
// Created by 梦想家xixi on 2022/1/4.
//

#ifndef MUTECT2CPP_MASTER_PEUTILS_H
#define MUTECT2CPP_MASTER_PEUTILS_H

#include "htslib/sam.h"
#include "cigar/Cigar.h"

class PeUtils {
public:
    bool isBeforeSoftClip();
    bool isImmediatelyBefore(CigarOperator cigarOperator);
    PeUtils(bam1_t* pe, int pos);
    bool isDeletion();
    bool isBeforeDeletionStart();
    bool isBeforeInsertion();
    bool isOnGenomeCigar(CigarOperator cigarOperator);
    int getLengthOfImmediatelyFollowingIndel();
    CigarElement & getCurrentCigarElement();
    CigarElement getNearestOnGenomeCigarElement(int direction);

private:
    int Cigar_offset;
    int pos;
    int currentStart;
    std::vector<CigarElement> nCigarElements;
    CigarElement currentCigarElement;
    bam1_t *pe;
    CigarElement getNextIndelCigarElement();
};


#endif //MUTECT2CPP_MASTER_PEUTILS_H
