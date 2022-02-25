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
    static std::shared_ptr<Cigar> calculateCigar(std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> altSeq, int altLength);
    static const SWParameters NEW_SW_PARAMETERS;
    static std::shared_ptr<Cigar> leftAlignCigarSequentially(std::shared_ptr<Cigar> & cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex);
    static bool isGood(const std::shared_ptr<Cigar> & c);
    static bool containsNOperator(std::vector<CigarElement> cigarElements);

private:
    static const int SW_PAD = 10;
    static bool isSWFailure(SmithWatermanAlignment* alignment);
    static bool hasConsecutiveIndels(std::vector<CigarElement> & elems);
    static bool startsWithDeletionIgnoringClips(std::vector<CigarElement> & elems);
};


#endif //MUTECT2CPP_MASTER_CIGARUTILS_H
