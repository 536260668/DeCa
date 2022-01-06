//
// Created by 梦想家xixi on 2022/1/4.
//

#ifndef MUTECT2CPP_MASTER_PEUTILS_H
#define MUTECT2CPP_MASTER_PEUTILS_H

#include "htslib/sam.h"
#include "cigar/Cigar.h"

class PeUtils {
public:
    static const char  DELETION_QUAL = 16;
    bool isBeforeSoftClip();
    bool isImmediatelyBefore(CigarOperator cigarOperator);
    bool isImmediatelyAfter(CigarOperator op);
    bool isAfterSoftClip();
    PeUtils(bam1_t* pe, int pos);
    bool isDeletion();
    bool isBeforeDeletionStart();
    bool isBeforeInsertion();
    bool isOnGenomeCigar(CigarOperator cigarOperator);
    int getLengthOfImmediatelyFollowingIndel();
    CigarElement & getCurrentCigarElement();
    CigarElement getNearestOnGenomeCigarElement(int direction);
    uint8_t getQual();
    uint8_t getBaseQuality(int pos);

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
