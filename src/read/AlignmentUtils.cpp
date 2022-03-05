//
// Created by 梦想家xixi on 2021/11/10.
//

#include "AlignmentUtils.h"

#include <utility>
#include <cstring>
#include "Mutect2Utils.h"

std::shared_ptr<Cigar> AlignmentUtils::consolidateCigar(std::shared_ptr<Cigar> c) {
    if(c == nullptr) {throw std::invalid_argument("Cigar cannot be null");}

    if(!needsConsolidation(c))
        return c;

    std::shared_ptr<Cigar> returnCigar(new Cigar());
    int sumLength = 0;
    CigarElement* lastElement = nullptr;

    for(CigarElement & cur : c->getCigarElements()) {
        if( cur.getLength() == 0)
            continue;

        if(lastElement != nullptr && lastElement->getOperator() != cur.getOperator()) {
            returnCigar->add(CigarElement(sumLength, lastElement->getOperator()));
            sumLength = 0;
        }
        sumLength += cur.getLength();
        lastElement = &cur;
    }
    if(sumLength > 0) {
        returnCigar->add(CigarElement(sumLength, lastElement->getOperator()));
    }
    return returnCigar;
}

bool AlignmentUtils::needsConsolidation(const std::shared_ptr<Cigar>& c) {
    if(c->numCigarElements() <= 1)
        return false;
    CigarOperator lastOp;
    for(CigarElement cur : c->getCigarElements()) {
        if(cur.getLength() == 0 || lastOp == cur.getOperator())
            return true;
        lastOp = cur.getOperator();
    }
    return false;
}

std::shared_ptr<uint8_t[]> AlignmentUtils::getBasesCoveringRefInterval(int refStart, int refEnd, std::shared_ptr<uint8_t[]>bases, int length, int basesStartOnRef,
                                                     const std::shared_ptr<Cigar>& basesToRefCigar) {
    if(refStart < 0 || refEnd < refStart) {
        throw std::invalid_argument("Bad start and/or stop");
    }
    Mutect2Utils::validateArg(basesStartOnRef >= 0, "BasesStartOnRef must be >= 0");
    Mutect2Utils::validateArg(bases != nullptr, "Null is not allowed there");
    Mutect2Utils::validateArg(basesToRefCigar != nullptr, "Null is not allowed there");
    Mutect2Utils::validateArg(length == basesToRefCigar->getReadLength(), "Mismatch in length between reference bases and cigar length.");

    int refPos = basesStartOnRef;
    int basesPos = 0;
    int basesStart = -1;
    int basesStop = -1;
    bool done = false;

    for(int iii = 0; ! done && iii <basesToRefCigar->numCigarElements(); iii++) {
        CigarElement ce = basesToRefCigar->getCigarElement(iii);
        switch (ce.getOperator()) {
            case I:
                basesPos += ce.getLength();
                break;
            case M:
            case X:
            case EQ:
                for (int i = 0; i < ce.getLength(); i++) {
                    if(refPos == refStart)
                        basesStart = basesPos;
                    if(refPos == refEnd) {
                        basesStop = basesPos;
                        done = true;
                        break;
                    }
                    refPos++;
                    basesPos++;
                }
                break;
            case D:
                for(int i = 0; i < ce.getLength(); i++) {
                    if(refPos == refEnd || refPos == refStart) {
                        return nullptr;
                    }
                    refPos++;
                }
                break;
            default:
                throw std::invalid_argument("Unsupported operator");
        }
    }

    if(basesStart == -1 || basesStop == -1)
        throw std::invalid_argument("Never found start or stop");

    std::shared_ptr<uint8_t[]> ret(new uint8_t[basesStop - basesStart + 2]);
    memcpy(ret.get(), bases.get() + basesStart, basesStop - basesStart + 2);
    return ret;
}

std::shared_ptr<Cigar> AlignmentUtils::trimCigarByReference(const std::shared_ptr<Cigar>&cigar, const int start, const int end) {
    Mutect2Utils::validateArg(start >= 0, "Start must be >= 0");
    Mutect2Utils::validateArg(end >= start, "End is < start");
    Mutect2Utils::validateArg(end <= cigar->getReferenceLength(), "End is beyond the cigar's reference length");

    std::shared_ptr<Cigar> result = trimCigar(cigar, start, end, true);
    Mutect2Utils::validateArg(result->getReferenceLength() == end - start + 1, "trimCigarByReference failure");
    return result;
}

std::shared_ptr<Cigar> AlignmentUtils::trimCigar(const std::shared_ptr<Cigar>& cigar, const int start, const int end, const bool byReference) {
    std::vector<CigarElement> newElements;

    int pos = 0;
    for (CigarElement elt : cigar->getCigarElements()) {
        if ( pos > end && (byReference || elt.getOperator() != D) ) break;

        switch (elt.getOperator()) {
            case D:
                if(!byReference) {
                    if(pos >= start)
                        newElements.emplace_back(elt);
                    break;
                }
            case EQ:
            case M:
            case X:
                pos = addCigarElements(newElements, pos, start, end, elt);
                break;
            case S:
            case I:
                if(byReference) {
                    if(pos >= start)
                        newElements.emplace_back(elt);
                } else {
                    pos = addCigarElements(newElements, pos, start, end, elt);
                }
                break;
            default:
                throw std::invalid_argument("Cannot handle");
        }
    }
    std::shared_ptr<Cigar> tmp(new Cigar(newElements));
    std::shared_ptr<Cigar> ret = consolidateCigar(tmp);
    return ret;
}

int AlignmentUtils::addCigarElements(std::vector<CigarElement> & dest, int pos, int start, int end, CigarElement elt) {
    int length = std::min(pos + elt.getLength() - 1, end) - std::max(pos, start) + 1;
    if(length > 0)
        dest.emplace_back(CigarElement(length, elt.getOperator()));
    return pos + elt.getLength();
}

bool AlignmentUtils::startsOrEndsWithInsertionOrDeletion(const std::shared_ptr<Cigar>& cigar) {
    Mutect2Utils::validateArg(cigar != nullptr, "Null is not allowed");
    if(cigar->isEmpty())
        return false;
    CigarOperator first = cigar->getCigarElement(0).getOperator();
    CigarOperator last = cigar->getCigarElement(cigar->numCigarElements()-1).getOperator();
    return first == D || first == I || last == D || last == I;
}

std::shared_ptr<Cigar> AlignmentUtils::removeTrailingDeletions(std::shared_ptr<Cigar> c) {
    std::vector<CigarElement> elements = c->getCigarElements();
    if(elements.at(elements.size()-1).getOperator() != D)
        return c;
    std::vector<CigarElement> newElements(elements.begin(), elements.end()-1);
    return std::shared_ptr<Cigar>(new Cigar(newElements));
}

std::shared_ptr<Cigar> AlignmentUtils::trimCigarByBases(const std::shared_ptr<Cigar>&cigar, int start, int end) {
    Mutect2Utils::validateArg(start >= 0, "tart must be >= 0 but got");
    Mutect2Utils::validateArg(end >= start, "End is < start ");
    Mutect2Utils::validateArg(end <= (cigar->getReadLength()), "End is beyond the cigar's read length");

    std::shared_ptr<Cigar> result = trimCigar(cigar, start, end, false);
    int expectedSize = end - start + 1;
    Mutect2Utils::validateArg((result->getReadLength()) == expectedSize, "trimCigarByBases failure");
    return result;
}

std::shared_ptr<Cigar>
AlignmentUtils::leftAlignSingleIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]>refSeq, int refLength, std::shared_ptr<uint8_t[]>readSeq, int readLength,
                                     int refIndex, int readIndex, bool cleanupCigar) {
    return leftAlignSingleIndel(std::move(cigar), refSeq, refLength, readSeq, readLength, refIndex, readIndex, 0, false);
}

std::shared_ptr<Cigar>
AlignmentUtils::leftAlignSingleIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]>refSeq, int refLength, std::shared_ptr<uint8_t[]>readSeq, int readLength,
                                     int refIndex, int readIndex, int leftmostAllowedAlignment, bool cleanupCigar1) {
    ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

    int indexOfIndel = -1;
    for(int i = 0; i < cigar->numCigarElements(); i++) {
        CigarElement ce = cigar->getCigarElement(i);
        if(ce.getOperator() == D || ce.getOperator() == I) {
            if(indexOfIndel != -1)
                throw std::invalid_argument("attempting to left align a CIGAR that has more than 1 indel in its alignment");
            indexOfIndel = i;
        }
    }

    if(indexOfIndel == -1)
        throw std::invalid_argument("attempting to left align a CIGAR that has no indels in its alignment");
    if(indexOfIndel == 0)
        return cigar;

    int indelLength = cigar->getCigarElement(indexOfIndel).getLength();
    int altStringLength = 0;
    std::shared_ptr<uint8_t[]> altString = createIndelString(cigar, indexOfIndel, refSeq, refLength, readSeq, readLength, refIndex, refIndex, altStringLength);
    if(altString == nullptr)
        return cigar;
    std::shared_ptr<Cigar> newCigar(new Cigar(*cigar));
    for(int i = 0; i < indelLength; i++) {
        std::shared_ptr<Cigar> tmp = moveCigarLeft(newCigar, indexOfIndel);
        newCigar = tmp;
        if(isIndelAlignedTooFarLeft(newCigar,leftmostAllowedAlignment)){
            break;
        }
        int newAltStringLength = 0;
        std::shared_ptr<uint8_t[]> newAltString = createIndelString(newCigar, indexOfIndel, refSeq, refLength, readSeq, readLength, refIndex, readIndex, newAltStringLength);
        bool reachedEndOfRead = cigarHasZeroSizeElement(newCigar);
        bool isEqual = true;
        if(altStringLength != newAltStringLength)
            isEqual = false;
        else {
            uint8_t * alt = altString.get();
            uint8_t * newAlt = newAltString.get();
            for(int j = 0; j < altStringLength; j++) {
                if(alt[j] != newAlt[j]){
                    isEqual = false;
                    break;
                }
            }
        }
        if(isEqual) {
            cigar = newCigar;
            i = -1;
            if(reachedEndOfRead)
                cigar = cleanupCigar1 ? cleanUpCigar(cigar) : cigar;
        }
        if(reachedEndOfRead)
            break;
    }
    return cigar;
}

void AlignmentUtils::ensureLeftAlignmentHasGoodArguments(const std::shared_ptr<Cigar>&cigar, std::shared_ptr<uint8_t[]>refSeq, std::shared_ptr<uint8_t[]>readSeq,
                                                          int refIndex, int readIndex) {
    Mutect2Utils::validateArg(cigar.get(), "ciagr");
    Mutect2Utils::validateArg(refSeq.get(), "refSeq");
    Mutect2Utils::validateArg(readSeq.get(), "readSeq");
    Mutect2Utils::validateArg(refIndex >= 0, "attempting to left align with a reference index less than 0");
    Mutect2Utils::validateArg(readIndex >= 0, "attempting to left align with a read index less than 0");
}

std::shared_ptr<uint8_t[]>
AlignmentUtils::createIndelString(const std::shared_ptr<Cigar>& cigar, int indexOfIndel, std::shared_ptr<uint8_t[]>refSeq, int refLength, std::shared_ptr<uint8_t[]>readSeq,
                                  int readLength, int refIndex, int readIndex, int & newLength) {
    CigarElement indel = cigar->getCigarElement(indexOfIndel);
    int indelLength = indel.getLength();

    int totalRefBases = 0;
    for(int i = 0; i < indexOfIndel; i++) {
        CigarElement ce = cigar->getCigarElement(i);
        int length = ce.getLength();

        switch (ce.getOperator()) {
            case M:
            case EQ:
            case X:
                readIndex += length;
                refIndex += length;
                totalRefBases += length;
                break;
            case S:
                readIndex += length;
                break;
            case N:
                refIndex += length;
                totalRefBases += length;
                break;
            default:
                break;
        }
    }

    if(totalRefBases + indelLength > refLength)
        indelLength -= (totalRefBases + indelLength - refLength);

    int altLength = refLength + (indelLength * (indel.getOperator() == D ? -1 : 1));
    std::shared_ptr<uint8_t[]> alt{new uint8_t[altLength]};

    if(refIndex > altLength || refIndex > refLength)
        return nullptr;
    memcpy(alt.get(), refSeq.get(), refIndex);
    int currentPos = refIndex;

    if(indel.getOperator() == D) {
        refIndex += indelLength;
    } else {
        memcpy(alt.get()+currentPos, readSeq.get()+readIndex, indelLength);
        currentPos += indelLength;
    }

    if(refLength - refIndex > altLength - currentPos)
        return nullptr;
    memcpy(alt.get()+currentPos, refSeq.get()+refIndex, refLength-refIndex);
    newLength = altLength;
    return alt;
}

std::shared_ptr<Cigar> AlignmentUtils::moveCigarLeft(const std::shared_ptr<Cigar>&cigar, int indexOfIndel) {
    std::vector<CigarElement> elements;
    elements.reserve(indexOfIndel - 1);
    for(int i = 0; i < indexOfIndel - 1; i++)
        elements.emplace_back(cigar->getCigarElement(i));

    CigarElement ce = cigar->getCigarElement(indexOfIndel - 1);
    elements.emplace_back(CigarElement(std::max(ce.getLength()-1, 0), ce.getOperator()));
    elements.emplace_back(cigar->getCigarElement(indexOfIndel));
    if(indexOfIndel + 1 < cigar->numCigarElements()) {
        ce = cigar->getCigarElement(indexOfIndel + 1);
        elements.emplace_back(CigarElement(ce.getLength() + 1, ce.getOperator()));
    } else {
        elements.emplace_back(CigarElement(1, M));
    }

    for (int i = indexOfIndel + 2; i < cigar->numCigarElements(); i++)
        elements.emplace_back(cigar->getCigarElement(i));

    return std::shared_ptr<Cigar>(new Cigar(elements));
}

bool AlignmentUtils::isIndelAlignedTooFarLeft(const std::shared_ptr<Cigar>&cigar, int leftmostAllowedAlignment) {
    int location = 0;
    for(CigarElement element : cigar->getCigarElements()) {
        if(element.getOperator() == D || element.getOperator() == I) {
            return location<leftmostAllowedAlignment;
        }
        if(CigarOperatorUtils::getConsumesReferenceBases(element.getOperator())) {
            location += element.getLength();
        }
    }
    return false;
}

bool AlignmentUtils::cigarHasZeroSizeElement(const std::shared_ptr<Cigar>&c) {
    for(CigarElement ce : c->getCigarElements()) {
        if(ce.getLength() == 0)
            return true;
    }
    return false;
}

std::shared_ptr<Cigar> AlignmentUtils::cleanUpCigar(const std::shared_ptr<Cigar>& c) {
    std::vector<CigarElement> elements;
    for(CigarElement ce : c->getCigarElements()) {
        if(ce.getLength() != 0 && (! elements.empty() || ce.getOperator() != D)) {
            elements.emplace_back(ce);
        }
    }

    return std::shared_ptr<Cigar>(new Cigar(elements));
}
