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
    static std::shared_ptr<Cigar> calculateCigar(const std::shared_ptr<uint8_t[]>& refSeq, int refLength, const std::shared_ptr<uint8_t[]>& altSeq, int altLength);
    static const SWParameters NEW_SW_PARAMETERS;
    static std::shared_ptr<Cigar> leftAlignCigarSequentially(std::shared_ptr<Cigar> & cigar, const std::shared_ptr<uint8_t[]>& refSeq, int refLength, const std::shared_ptr<uint8_t[]>& readSeq, int readLength, int refIndex, int readIndex);
    static bool isGood(const std::shared_ptr<Cigar> & c);
    static bool isGood(bam1_t * read);
    static bool containsNOperator(std::vector<CigarElement> cigarElements);
    static bool containsNOperator(uint32_t n_cigar, uint32_t* cigarArray);

private:
    static const int SW_PAD = 10;
    static bool isSWFailure(SmithWatermanAlignment* alignment);
    static bool hasConsecutiveIndels(std::vector<CigarElement> & elems);
    static bool hasConsecutiveIndels(uint32_t * cigarArray, uint32_t n_cigar);
    static bool startsWithDeletionIgnoringClips(const std::vector<CigarElement> & elems);
    static bool startsWithDeletionIgnoringClips(uint32_t * cigarArray, uint32_t n_cigar);
};


#endif //MUTECT2CPP_MASTER_CIGARUTILS_H
