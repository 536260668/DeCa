//
// Created by 梦想家xixi on 2021/12/21.
//

#include "ClippingOp.h"
#include <stack>
#include "ReadUtils.h"

ClippingOp::ClippingOp(int start, int stop) : start(start), stop(stop){
}

std::shared_ptr<SAMRecord> ClippingOp::applyHardClipBases(std::shared_ptr<SAMRecord> read, int start, int stop) {
    Cigar* cigar = read->getCigar();
    Cigar* tmp = nullptr;
    CigarShift* cigarShift;
    if(read->isUnmapped()) {
        tmp = new Cigar();
        cigarShift = new CigarShift(tmp, 0, 0);
    } else {
        cigarShift = hardClipCigar(cigar, start, stop);
    }
    int newLength = read->getLength() - (stop - start + 1) - cigarShift->shiftFromStart - cigarShift->shiftFromEnd;
    if(newLength == 0) {
        delete cigarShift->cigar;
        delete cigarShift;
        return ReadUtils::emptyRead(read);
    }
    uint8_t * newBases = new uint8_t[newLength];
    uint8_t * newQuals = new uint8_t[newLength];
    int copyStart = (start == 0) ? stop + 1 + cigarShift->shiftFromStart : cigarShift->shiftFromStart;
    memcpy(newBases, read->getBasesNoCopy()+copyStart, newLength);
    memcpy(newQuals, read->getBaseQualitiesNoCopy()+copyStart, newLength);
    std::shared_ptr<SAMRecord> hardClippedRead{new SAMRecord(*read)};
    hardClippedRead->setBaseQualities(newQuals, newLength);
    hardClippedRead->setBases(newBases, newLength);
    hardClippedRead->setCigar(cigarShift->cigar);
    if(start == 0 && !read->isUnmapped()) {
        hardClippedRead->setPosition(read->getContig(), read->getStart() + calculateAlignmentStartShift(cigar, cigarShift->cigar));
    }
    if(ReadUtils::hasBaseIndelQualities(read)) {
        uint8_t * newBaseInsertionQuals = new uint8_t[newLength];
        uint8_t * newBaseDeletionQuals = new uint8_t[newLength];
        int length;
        uint8_t * new_base = ReadUtils::getBaseInsertionQualities(read, length);
        memcpy(newBaseInsertionQuals, new_base+copyStart, newLength);
        delete[] new_base;
        new_base = ReadUtils::getBaseDeletionQualities(read, length);
        memcpy(newBaseDeletionQuals, new_base+copyStart, newLength);
        delete[] new_base;
        ReadUtils::setDeletionBaseQualities(hardClippedRead, newBaseDeletionQuals, newLength);
        ReadUtils::setInsertionBaseQualities(hardClippedRead, newBaseInsertionQuals, newLength);
    }
    return hardClippedRead;
}

ClippingOp::CigarShift *ClippingOp::hardClipCigar(Cigar *cigar, int start, int stop) {
    Cigar* newCigar = new Cigar();
    int index = 0;
    int totalHardClipCount = stop - start + 1;
    int alignmentShift = 0;

    if(start == 0) {
        std::vector<CigarElement> cigarElementIterator = cigar->getCigarElements();
        int i = 0;
        CigarElement cigarElement = cigarElementIterator[0];
        i++;
        while(cigarElement.getOperator() == H) {
            totalHardClipCount += cigarElement.getLength();
            if(i < cigarElementIterator.size()) {
                cigarElement = cigarElementIterator[i];
                i++;
            } else {
                throw std::invalid_argument("Read is entirely hard-clipped, shouldn't be trying to clip it's cigar string");
            }
        }
        while(index <= stop) {
            int shift = 0;
            if(CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                shift = cigarElement.getLength();
            }
            if(index + shift == stop + 1) {
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());
                newCigar->add(CigarElement(totalHardClipCount + alignmentShift, H));
            } else if (index + shift > stop + 1) {
                int elementLengthAfterChopping = cigarElement.getLength() - (stop - index + 1);
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, stop - index + 1);
                newCigar->add(CigarElement(totalHardClipCount + alignmentShift, H));
                newCigar->add(CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
            }
            index += shift;
            alignmentShift += calculateHardClippingAlignmentShift(cigarElement, shift);

            if(index <= stop && i < cigarElementIterator.size()) {
                cigarElement = cigarElementIterator[i];
                i++;
            } else {
                break;
            }
        }

        while(i < cigarElementIterator.size()) {
            cigarElement = cigarElementIterator[i];
            newCigar->add(CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
        }
    } else {
        std::vector<CigarElement> cigarElementIterator = cigar->getCigarElements();
        int i = 0;
        CigarElement cigarElement = cigarElementIterator[i];
        i++;

        while(index < start) {
            int shift = 0;
            if(CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                shift = cigarElement.getLength();
            }

            if(index + shift < start) {
                newCigar->add(CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
            } else {
                int elementLengthAfterChopping = start - index;
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength() - (start - index));

                if(cigarElement.getOperator() == H) {
                    totalHardClipCount += elementLengthAfterChopping;
                } else {
                    newCigar->add(CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
            }
            index += shift;
            if(index < start && i < cigarElementIterator.size()) {
                cigarElement = cigarElementIterator[i];
                i++;
            } else {
                break;
            }
        }
        while(i < cigarElementIterator.size()) {
            cigarElement = cigarElementIterator[i];
            i++;
            alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());

            if(cigarElement.getOperator() == H) {
                totalHardClipCount += cigarElement.getLength();
            }
        }
        newCigar->add(CigarElement(totalHardClipCount + alignmentShift, H));
    }
    return cleanHardClippedCigar(newCigar);
}

int ClippingOp::calculateHardClippingAlignmentShift(CigarElement &cigarElement, int clippedLength) {
    if(cigarElement.getOperator() == I) {
        return -clippedLength;
    } else if (cigarElement.getOperator() == D || cigarElement.getOperator() == S) {
        return cigarElement.getLength();
    }
    return 0;
}

ClippingOp::CigarShift *ClippingOp::cleanHardClippedCigar(Cigar *cigar) {
    Cigar* cleanCigar = new Cigar();
    int shiftFromStart = 0;
    int shiftFromEnd = 0;
    std::stack<CigarElement> cigarStack;
    std::stack<CigarElement> inverseCigarStack;
    for(CigarElement cigarElement : cigar->getCigarElements()) {
        cigarStack.push(cigarElement);
    }
    for(int i = 1; i <= 2; i++) {
        int shift = 0;
        int totalHardClip = 0;
        bool readHasStarted = false;
        bool addedHardClips = false;

        while(!cigarStack.empty()) {
            CigarElement cigarElement = cigarStack.top();
            cigarStack.pop();
            if(!readHasStarted &&
            cigarElement.getOperator() != D &&
            cigarElement.getOperator() != N &&
            cigarElement.getOperator() != H) {
                readHasStarted = true;
            } else if (!readHasStarted && cigarElement.getOperator() == H) {
                totalHardClip += cigarElement.getLength();
            } else if (!readHasStarted && cigarElement.getOperator() == D) {
                totalHardClip += cigarElement.getLength();
            } else if (!readHasStarted && cigarElement.getOperator() == N) {
                totalHardClip += cigarElement.getLength();
            }

            if(readHasStarted) {
                if(i == 1) {
                    if(!addedHardClips) {
                        if(totalHardClip > 0) {
                            inverseCigarStack.push(CigarElement(totalHardClip, H));
                        }
                        addedHardClips = true;
                    }
                    inverseCigarStack.push(cigarElement);
                } else {
                    if(!addedHardClips) {
                        if(totalHardClip > 0) {
                            cleanCigar->add(CigarElement(totalHardClip, H));
                        }
                        addedHardClips = true;
                    }
                    cleanCigar->add(cigarElement);
                }
            }
        }
        if(i == 1) {
            shiftFromEnd = shift;
            cigarStack = inverseCigarStack;
        } else {
            shiftFromStart = shift;
        }
    }
    delete cigar;
    return new CigarShift(cleanCigar, shiftFromStart, shiftFromEnd);
}

int ClippingOp::calculateAlignmentStartShift(Cigar *oldCigar, Cigar *newCigar) {
    int newShift = calcHardSoftOffset(newCigar);
    int oldShift = calcHardSoftOffset(oldCigar);
    return newShift - oldShift;
}

int ClippingOp::calcHardSoftOffset(Cigar *cigar) {
    std::vector<CigarElement> elements = cigar->getCigarElements();
    int size = 0;
    int i = 0;
    while ( i < elements.size() && elements[i].getOperator() == H ) {
        size += elements[i].getLength();
        i++;
    }
    while ( i < elements.size() && elements[i].getOperator() == S ) {
        size += elements[i].getLength();
        i++;
    }
    return size;
}

std::shared_ptr<SAMRecord> ClippingOp::apply(ClippingRepresentation algorithm, std::shared_ptr<SAMRecord> originalRead, bool runAsserts) {
    switch(algorithm){
        case HARDCLIP_BASES: {
            return applyHardClipBases(originalRead, start, stop);
        }
        default: {
            throw std::invalid_argument("Unexpected Clipping operator type");
        }
    }
}

ClippingOp::CigarShift::CigarShift(Cigar *cigar, int shiftFromStart, int shiftFromEnd) : cigar(cigar), shiftFromStart(shiftFromStart), shiftFromEnd(shiftFromEnd){
}
