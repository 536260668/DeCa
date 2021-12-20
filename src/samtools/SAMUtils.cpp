//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMUtils.h"

int SAMUtils::getUnclippedStart(int alignmentStart, Cigar *cigar) {
    int unClippedStart = alignmentStart;
    for(CigarElement cig : cigar->getCigarElements()) {
        CigarOperator op = cig.getOperator();
        if (op != S && op != H) {
            break;
        }
        unClippedStart -= cig.getLength();
    }
    return unClippedStart;
}

int SAMUtils::getUnclippedEnd(int alignmentEnd, Cigar *cigar) {
    int unClippedEnd = alignmentEnd;
    std::vector<CigarElement> cigs = cigar->getCigarElements();
    for(int i = cigs.size() - 1; i >= 0; --i) {
        CigarElement cig = cigs[i];
        CigarOperator op = cig.getOperator();
        if (op != S && op != H) {
            break;
        }

        unClippedEnd += cig.getLength();
    }

    return unClippedEnd;

}

bool SAMUtils::isValidUnsignedIntegerAttribute(long value) {
    return value >= 0L && value <= 4294967295L;
}

short SAMUtils::makeBinaryTag(std::string &tag) {
    if(tag.length() != 2) {
        throw std::invalid_argument("String tag does not have length() == 2: ");
    }
    char a = tag[1];
    char b = tag[0];
    return a << 8 | b;
}
