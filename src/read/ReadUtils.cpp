//
// Created by 梦想家xixi on 2021/12/18.
//

#include "ReadUtils.h"
#include "samtools/SAMUtils.h"

std::string ReadUtils::BQSR_BASE_DELETION_QUALITIES = "BD";
std::string ReadUtils::BQSR_BASE_INSERTION_QUALITIES = "BI";

SAMRecord *ReadUtils::emptyRead(SAMRecord *read) {
    SAMRecord* emptyRead = new SAMRecord(*read);
    emptyRead->setIsUnmapped();
    emptyRead->setMappingQuality(0);
    emptyRead->setCigar(new Cigar());
    emptyRead->setBases(nullptr, 0);
    emptyRead->setBaseQualities(nullptr, 0);
    emptyRead->clearAttributes();

    std::string readGroup = read->getReadGroup();
    if(readGroup.empty()) {
        emptyRead->setAttribute((std::string&)"RG", readGroup);
    }
    return emptyRead;
}

void ReadUtils::assertAttributeNameIsLegal(std::string &attributeName) {
    if(attributeName.empty() || attributeName.length() !=2) {
        throw std::invalid_argument("Read attribute invalid: attribute names must be non-null two-character Strings matching the pattern /[A-Za-z][A-Za-z0-9]/");
    }
}

int
ReadUtils::getReadCoordinateForReferenceCoordinate(int alignmentStart, Cigar *cigar, int refCoord, ClippingTail tail,
                                                   bool allowGoalNotReached) {
    std::pair<int, bool> result = getReadCoordinateForReferenceCoordinate(alignmentStart, cigar, refCoord, allowGoalNotReached);
    int readCoord = result.first;
    if(result.second && tail == RIGHT_TAIL){
        readCoord++;
    }
    CigarElement* firstElementIsInsertion = readStartsWithInsertion(cigar);
    if(readCoord == 0 && tail == LEFT_TAIL && firstElementIsInsertion != nullptr) {
        readCoord = std::min(firstElementIsInsertion->getLength(), cigar->getReadLength() - 1);
    }
    return readCoord;
}

std::pair<int, bool> ReadUtils::getReadCoordinateForReferenceCoordinate(int alignmentStart, Cigar *cigar, int refCoord,
                                                                        bool allowGoalNotReached) {
    int readBases = 0;
    int refBases = 0;
    bool fallsInsideDeletionOrSkippedRegion = false;
    bool endJustBeforeDeletionOrSkippedRegion = false;
    bool fallsInsideOrJustBeforeDeletionOrSkippedRegion = false;
    int goal = refCoord - alignmentStart;

    if (goal < 0) {
        if (allowGoalNotReached) {
            return {CLIPPING_GOAL_NOT_REACHED, false};
        } else {
            throw std::invalid_argument("Somehow the requested coordinate is not covered by the read. Too many deletions?");
        }
    }
    bool goalReached = refBases == goal;
    std::vector<CigarElement> cigars = cigar->getCigarElements();
    for(int i = 0; !goalReached && i < cigars.size();) {
        CigarElement cigarElement = cigars[i];
        i++;
        int shift = 0;
        if(CigarOperatorUtils::getConsumesReferenceBases(cigarElement.getOperator()) || cigarElement.getOperator() == S) {
            if (refBases + cigarElement.getLength() < goal) {
                shift = cigarElement.getLength();
            } else {
                shift = goal - refBases;
            }

            refBases += shift;
        }
        goalReached = refBases == goal;
        if(!goalReached && CigarOperatorUtils::getConsumesReferenceBases(cigarElement.getOperator())) {
            refBases += cigarElement.getLength();
        }
        if(goalReached) {
            bool endsWithinCigar = shift < cigarElement.getLength();
            if(!endsWithinCigar && i >= cigars.size()) {
                if(allowGoalNotReached) {
                    return {CLIPPING_GOAL_NOT_REACHED, false};
                } else {
                    throw std::invalid_argument("Reference coordinate corresponds to a non-existent base in the read.");
                }
            }
            CigarElement* nextCigarElement = nullptr;
            if(endsWithinCigar) {
                fallsInsideDeletionOrSkippedRegion = (cigarElement.getOperator() == D || cigarElement.getOperator() == N);
            } else {
                nextCigarElement = &cigars[i];
                i++;
                if(nextCigarElement->getOperator() == I) {
                    readBases += nextCigarElement->getLength();
                    if(i >= cigars.size()) {
                        if(allowGoalNotReached) {
                            return {CLIPPING_GOAL_NOT_REACHED, false};
                        } else {
                            throw std::invalid_argument("Reference coordinate corresponds to a non-existent base in the read.");
                        }
                    }
                    nextCigarElement = &cigars[i];
                    i++;
                }
                endJustBeforeDeletionOrSkippedRegion = (nextCigarElement->getOperator() == D || nextCigarElement->getOperator() == N);
            }
            fallsInsideOrJustBeforeDeletionOrSkippedRegion = endJustBeforeDeletionOrSkippedRegion || fallsInsideDeletionOrSkippedRegion;
            if(!fallsInsideOrJustBeforeDeletionOrSkippedRegion && CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                readBases += shift;
            } else if(endJustBeforeDeletionOrSkippedRegion && CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                readBases += shift - 1;
            } else if(fallsInsideDeletionOrSkippedRegion ||
                    (endJustBeforeDeletionOrSkippedRegion && nextCigarElement->getOperator() == N) ||
                    (endJustBeforeDeletionOrSkippedRegion && nextCigarElement->getOperator() == D)) {
                readBases--;
            }
        }
    }
    if(!goalReached) {
        if(allowGoalNotReached) {
            return {CLIPPING_GOAL_NOT_REACHED, false};
        } else {
            throw std::invalid_argument("Somehow the requested coordinate is not covered by the read. Alignment");
        }
    }

    return {readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion};
}

CigarElement* ReadUtils::readStartsWithInsertion(Cigar *cigarForRead, bool ignoreSoftClipOps) {
    for(CigarElement cigarElement : cigarForRead->getCigarElements()) {
        if(cigarElement.getOperator() == I) {
            return &cigarElement;
        } else if (cigarElement.getOperator() != H && (!ignoreSoftClipOps || cigarElement.getOperator() != S)) {
            break;
        }
    }
    return nullptr;
}

CigarElement* ReadUtils::readStartsWithInsertion(Cigar *cigarForRead) {
    return readStartsWithInsertion(cigarForRead, true);
}

int ReadUtils::getReadCoordinateForReferenceCoordinate(SAMRecord *read, int refCoord, ClippingTail tail) {
    int leftmostSafeVariantPosition = std::max(read->getSoftStart(), refCoord);
    return getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), leftmostSafeVariantPosition, tail,
                                                   false);
}

int ReadUtils::getSoftStart(SAMRecord *read) {
    Mutect2Utils::validateArg(read, "read");

    int softStart = read->getStart();
    for(CigarElement cig : read->getCigarElements()) {
        CigarOperator op = cig.getOperator();
        if(op == S) {
            softStart -= cig.getLength();
        } else if (op != H) {
            break;
        }
    }
    return softStart;
}

int ReadUtils::getSoftEnd(SAMRecord *read) {
    Mutect2Utils::validateArg(read, "read");

    bool foundAlignedBase = false;
    int softEnd = read->getEnd();
    std::vector<CigarElement> cigs = read->getCigarElements();
    for (int i = cigs.size() - 1; i >= 0; --i) {
        CigarElement cig = cigs[i];
        CigarOperator op = cig.getOperator();

        if (op == S){
            softEnd += cig.getLength();
        } else if (op != H) {
            foundAlignedBase = true;
            break;
        }
    }
    if( !foundAlignedBase ) {
        softEnd = read->getEnd();
    }
    return softEnd;
}

bool ReadUtils::hasBaseIndelQualities(SAMRecord *read) {
    return read->getAttribute(SAMUtils::makeBinaryTag(BQSR_BASE_INSERTION_QUALITIES)) || read->getAttribute(SAMUtils::makeBinaryTag(BQSR_BASE_DELETION_QUALITIES));
}

uint8_t *ReadUtils::getExistingBaseInsertionQualities(SAMRecord *read, int &length) {
    std::string str = read->getAttributeAsString(BQSR_BASE_INSERTION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

uint8_t *ReadUtils::getExistingBaseDeletionQualities(SAMRecord *read, int &length) {
    std::string str = read->getAttributeAsString(BQSR_BASE_DELETION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

uint8_t *ReadUtils::getBaseInsertionQualities(SAMRecord *read, int &length) {
    uint8_t * quals = getExistingBaseInsertionQualities(read, length);
    if(quals == nullptr) {
        length = read->getBaseQualitiesLength();
        quals = new uint8_t[length]{45};
    }
    return quals;
}

uint8_t *ReadUtils::getBaseDeletionQualities(SAMRecord *read, int &length) {
    uint8_t * quals = getExistingBaseDeletionQualities(read, length);
    if(quals == nullptr) {
        length = read->getBaseQualitiesLength();
        quals = new uint8_t[length]{45};
    }
    return quals;
}

void ReadUtils::setInsertionBaseQualities(SAMRecord *read, uint8_t *quals, int length) {
    read->setAttribute(BQSR_BASE_INSERTION_QUALITIES, quals == nullptr ? nullptr : SAMUtils::p)
}
