//
// Created by 梦想家xixi on 2021/11/25.
//

#include "CigarUtils.h"
#include "Mutect2Utils.h"
#include "SWNativeAlignerWrapper.h"
#include "AlignmentUtils.h"

const SWParameters  CigarUtils::NEW_SW_PARAMETERS = SWParameters(200, -150, -260, -11);

Cigar *CigarUtils::calculateCigar(uint8_t *refSeq, int refLength, uint8_t *altSeq, int altLength) {
    Mutect2Utils::validateArg(refSeq, "refSeq");
    Mutect2Utils::validateArg(altSeq, "altSeq");
    if(altLength == 0) {
        std::vector<CigarElement> elements;
        elements.emplace_back(CigarElement(refLength, D));
        return new Cigar(elements);
    }
    if(altLength == refLength) {
        int mismatchCount = 0;
        for (int n = 0; n < refLength && mismatchCount <= 2; n++) {
            mismatchCount += (altSeq[n] == refSeq[n] ? 0 : 1);
        }
        if(mismatchCount <= 2) {
            Cigar* matching = new Cigar();
            matching->add(CigarElement(refLength, M));
            return matching;
        }
    }
    Cigar* nonStandard;
    int paddedRefLength = refLength + 2*SW_PAD;
    uint8_t * paddedRef = new uint8_t[paddedRefLength];
    int paddedPathLength = altLength + 2*SW_PAD;
    uint8_t * paddedPath = new uint8_t[paddedPathLength];
    uint8_t tmp[10] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};
    memcpy(paddedRef, tmp, 10);
    memcpy(paddedRef+10, refSeq, refLength);
    memcpy(paddedRef+10+refLength, tmp, 10);
    memcpy(paddedPath, tmp, 10);
    memcpy(paddedPath+10, altSeq, altLength);
    memcpy(paddedPath+10+altLength, tmp, 10);
    SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
    SmithWatermanAlignment* alignment = wrapper.align(paddedRef, paddedRefLength, paddedPath, paddedPathLength,
                                                      const_cast<SWParameters *>(&NEW_SW_PARAMETERS), SOFTCLIP);
    if ( isSWFailure(alignment) ) {
        return nullptr;
    }

    int baseStart = SW_PAD;
    int baseEnd = paddedPathLength - SW_PAD - 1;
    nonStandard = AlignmentUtils::trimCigarByBases(alignment->getCigar(), baseStart, baseEnd);
    if(nonStandard->getReferenceLength() != refLength) {
        nonStandard->add(CigarElement(refLength = nonStandard->getReferenceLength(), D));
    }

    return leftAlignCigarSequentially(nonStandard, refSeq, refLength, altSeq, altLength, 0, 0);
}

bool CigarUtils::isSWFailure(SmithWatermanAlignment *alignment) {
    if(alignment->getAlignmentOffset() > 0) {
        return true;
    }
    for(CigarElement ce : alignment->getCigar()->getCigarElements()) {
        if(ce.getOperator() == S)
            return true;
    }
    return false;
}

Cigar *
CigarUtils::leftAlignCigarSequentially(Cigar *cigar, uint8_t *refSeq, int refLength, uint8_t *readSeq, int readLength,
                                       int refIndex, int readIndex) {
    Mutect2Utils::validateArg(cigar, "cigar null");
    Mutect2Utils::validateArg(refSeq, "refSeq null");
    Mutect2Utils::validateArg(readSeq, "readSeq null");

    Cigar* cigarToReturn = new Cigar();
    Cigar* cigarToAlign = new Cigar();

    for(int i = 0; i < cigar->numCigarElements(); i++) {
        CigarElement ce = cigar->getCigarElement(i);
        if(ce.getOperator() == D || ce.getOperator() == I) {
            cigarToAlign->add(ce);
            Cigar* leftAligned = AlignmentUtils::leftAlignSingleIndel(cigarToAlign, refSeq, refLength, readSeq, readLength, refIndex, readIndex,
                                                                      false);
            for(CigarElement toAdd : leftAligned->getCigarElements()) {cigarToReturn->add(toAdd);}
            refIndex += cigarToAlign->getReferenceLength();
            readIndex += cigarToAlign->getReadLength();
            delete cigarToAlign;
            cigarToAlign = new Cigar();
        } else {
            cigarToAlign->add(ce);
        }
    }
    if(!cigarToAlign->isEmpty()) {
        for(CigarElement toAdd : cigarToAlign->getCigarElements()) {
            cigarToReturn->add(toAdd);
        }
    }

    Cigar* result = AlignmentUtils::consolidateCigar(cigarToReturn);
    Mutect2Utils::validateArg(result->getReferenceLength() == cigar->getReferenceLength(), "leftAlignCigarSequentially failed to produce a valid CIGAR.");
    return result;
}
