//
// Created by 梦想家xixi on 2021/12/18.
//

#include "ReadUtils.h"
#include "samtools/SAMUtils.h"
#include <cmath>

std::string ReadUtils::BQSR_BASE_DELETION_QUALITIES = "BD";
std::string ReadUtils::BQSR_BASE_INSERTION_QUALITIES = "BI";
const int ReadUtils::CLIPPING_GOAL_NOT_REACHED = -1;

std::shared_ptr<SAMRecord> ReadUtils::emptyRead(std::shared_ptr<SAMRecord> & read) {
    std::shared_ptr<SAMRecord> emptyRead(new SAMRecord(*read));
    emptyRead->setIsUnmapped();
    emptyRead->setMappingQuality(0);
    emptyRead->setCigar(new Cigar());
    emptyRead->setBases(nullptr, 0);
    emptyRead->setBaseQualities(nullptr, 0);
    emptyRead->clearAttributes();

//    std::string readGroup = read->getReadGroup();
//    if(readGroup.empty()) {
//        emptyRead->setAttribute((std::string&)"RG", readGroup);
//    }
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
        if(!goalReached && CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
            readBases += cigarElement.getLength();
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

int ReadUtils::getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> & read, int refCoord, ClippingTail tail) {
    //int leftmostSafeVariantPosition = std::max(read->getSoftStart(), refCoord);
    return getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), refCoord, tail,
                                                   false);
}

int ReadUtils::getSoftStart(std::shared_ptr<SAMRecord> & read) {
    Mutect2Utils::validateArg(read != nullptr, "read");

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

int ReadUtils::getSoftEnd(std::shared_ptr<SAMRecord> & read) {
    Mutect2Utils::validateArg(read != nullptr, "read");

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

bool ReadUtils::hasBaseIndelQualities(std::shared_ptr<SAMRecord> & read) {
    return read->getAttribute(SAMUtils::makeBinaryTag(BQSR_BASE_INSERTION_QUALITIES)) || read->getAttribute(SAMUtils::makeBinaryTag(BQSR_BASE_DELETION_QUALITIES));
}

uint8_t *ReadUtils::getExistingBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    std::string str = read->getAttributeAsString(BQSR_BASE_INSERTION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

uint8_t *ReadUtils::getExistingBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    std::string str = read->getAttributeAsString(BQSR_BASE_DELETION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

uint8_t *ReadUtils::getBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    uint8_t * quals = getExistingBaseInsertionQualities(read, length);
    if(quals == nullptr) {
        length = read->getBaseQualitiesLength();
        quals = new uint8_t[length]{45};
    }
    return quals;
}

uint8_t *ReadUtils::getBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    uint8_t * quals = getExistingBaseDeletionQualities(read, length);
    if(quals == nullptr) {
        length = read->getBaseQualitiesLength();
        quals = new uint8_t[length]{45};
    }
    return quals;
}

void ReadUtils::setInsertionBaseQualities(std::shared_ptr<SAMRecord> & read, uint8_t *quals, int length) {
    read->setAttribute(BQSR_BASE_INSERTION_QUALITIES, quals == nullptr ? "" : SAMUtils::phredToFastq(quals, length));
}

void ReadUtils::setDeletionBaseQualities(std::shared_ptr<SAMRecord> & read, uint8_t *quals, int length) {
    read->setAttribute(BQSR_BASE_DELETION_QUALITIES, quals  == nullptr ? "" : SAMUtils::phredToFastq(quals, length));
}

int ReadUtils::getAssignedReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader *header) {
    return header->getSequenceIndex(read->getAssignedContig());
}

int ReadUtils::getSAMFlagsForRead(std::shared_ptr<SAMRecord> & read) {
    int samFlags = 0;

    if(read->isPaired()) {
        samFlags |= 1;
    }
    if(read->isProperlyPaired()) {
        samFlags |= 2;
    }
    if(read->isUnmapped()) {
        samFlags |= 4;
    }
    if(read->isPaired() && read->mateIsUnmapped()) {
        samFlags |= 8;
    }
    if(!read->isUnmapped() && read->isReverseStrand()) {
        samFlags |= 16;
    }
    if(read->isPaired() && ! read->mateIsUnmapped() && read->mateIsReverseStrand()) {
        samFlags |= 32;
    }
    if(read->isFirstOfPair()) {
        samFlags |= 64;
    }
    if(read->isSecondOfPair()) {
        samFlags |= 128;
    }
    if(read->isSecondaryAlignment()) {
        samFlags |= 256;
    }
    if(read->failsVendorQualityCheck()) {
        samFlags |= 512;
    }
    if(read->isDuplicate()) {
        samFlags |= 1024;
    }
    if(read->isSupplementaryAlignment()) {
        samFlags |= 2048;
    }
    return samFlags;
}

int ReadUtils::getMateReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader *header) {
    if(read->mateIsUnmapped()) {
        return -1;
    }
    return header->getSequenceIndex(read->getMateContig());
}

bool ReadUtils::alignmentAgreesWithHeader(SAMFileHeader *header, std::shared_ptr<SAMRecord> & read) {
    int referenceIndex = getReferenceIndex(read, header);

    if(! read->isUnmapped() && referenceIndex == SAMRecord::NO_ALIGNMENT_REFERENCE_INDEX) {
        return false;
    }
    SAMSequenceRecord contigHeader = header->getSequenceDictionary().getSequences()[referenceIndex];
    return read->isUnmapped() || read->getStart() <= contigHeader.getSequenceLength();
}

int ReadUtils::getReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader *header) {
    if(read->isUnmapped()) {
        return SAMRecord::NO_ALIGNMENT_REFERENCE_INDEX;
    }
    return header->getSequenceIndex(read->getContig());
}

bool ReadUtils::hasWellDefinedFragmentSize(std::shared_ptr<SAMRecord> & read) {
    if(read->getFragmentLength() == 0) {
        return false;
    }
    if( ! read->isPaired()) {
        return false;
    }
    if(read->isUnmapped() || read->mateIsUnmapped()) {
        return false;
    }
    if(read->isReverseStrand() == read->mateIsReverseStrand()) {
        return false;
    }
    if(read->isReverseStrand()) {
        return read->getEnd() > read->getMateStart();
    } else {
        return read->getStart() <= read->getMateStart() + read->getFragmentLength();
    }
}

int ReadUtils::getAdaptorBoundary(std::shared_ptr<SAMRecord> & read) {
    if(!hasWellDefinedFragmentSize(read)) {
        return INT32_MIN;
    } else if (read->isReverseStrand()) {
        return read->getMateStart() - 1;
    } else {
        int insertSize = std::abs(read->getFragmentLength());
        return read->getStart() + insertSize;
    }
}

bool ReadUtils::isBaseInsideAdaptor(std::shared_ptr<SAMRecord> & read, long basePos) {
    int adaptorBoundary = read->getAdaptorBoundary();
    if(adaptorBoundary == INT32_MIN || read->getFragmentLength() > 100)
        return false;
    return read->isReverseStrand() ? basePos <= adaptorBoundary : basePos >= adaptorBoundary;
}

bool ReadUtils::isInsideRead(std::shared_ptr<SAMRecord> &read, int referenceCoordinate) {
    return referenceCoordinate >= read->getStart() && referenceCoordinate <= read->getEnd();
}

