//
// Created by 梦想家xixi on 2021/11/10.
//

#include "AlignmentUtils.h"
#include "Mutect2Utils.h"

Cigar *AlignmentUtils::consolidateCigar(Cigar *c) {
    if(c == nullptr) {throw std::invalid_argument("Cigar cannot be null");}

    if(!needsConsolidation(c))
        return c;

    Cigar* returnCigar = new Cigar();
    int sumLength = 0;
    CigarElement* lastElement = nullptr;

    for(CigarElement cur : c->getCigarElements()) {
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

bool AlignmentUtils::needsConsolidation(Cigar *c) {
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

uint8_t *AlignmentUtils::getBasesCoveringRefInterval(int refStart, int refEnd, uint8_t *bases, int length, int basesStartOnRef,
                                                     Cigar *basesToRefCigar) {
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

    uint8_t * ret = new uint8_t[basesStop - basesStart + 2];
    memcpy(ret, bases + basesStart, basesStop - basesStart + 2);
    return ret;
}

Cigar *AlignmentUtils::trimCigarByReference(Cigar *cigar, const int start, const int end) {
    Mutect2Utils::validateArg(start >= 0, "Start must be >= 0");
    Mutect2Utils::validateArg(end >= start, "End is < start");
    Mutect2Utils::validateArg(end <= cigar->getReferenceLength(), "End is beyond the cigar's reference length");

    Cigar * result = trimCigar(cigar, start, end, true);
    Mutect2Utils::validateArg(result->getReferenceLength() == end - start + 1, "trimCigarByReference failure");
    return result;
}

Cigar *AlignmentUtils::trimCigar(Cigar *cigar, const int start, const int end, const bool byReference) {
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
    Cigar* tmp = new Cigar(newElements);
    Cigar* ret = consolidateCigar(tmp);
    if(tmp != ret)
        delete tmp;
    return ret;
}

int AlignmentUtils::addCigarElements(std::vector<CigarElement> & dest, int pos, int start, int end, CigarElement elt) {
    int length = std::min(pos + elt.getLength() - 1, end) - std::max(pos, start) + 1;
    if(length > 0)
        dest.emplace_back(CigarElement(length, elt.getOperator()));
    return pos + elt.getLength();
}

bool AlignmentUtils::startsOrEndsWithInsertionOrDeletion(Cigar *cigar) {
    Mutect2Utils::validateArg(cigar != nullptr, "Null is not allowed");
    if(cigar->isEmpty())
        return false;
    CigarOperator first = cigar->getCigarElement(0).getOperator();
    CigarOperator last = cigar->getCigarElement(cigar->numCigarElements()-1).getOperator();
    return first == D || first == I || last == D || last == I;
}

Cigar *AlignmentUtils::removeTrailingDeletions(Cigar *c) {
    std::vector<CigarElement> elements = c->getCigarElements();
    if(elements.at(elements.size()-1).getOperator() != D)
        return c;
    std::vector<CigarElement> newElements(elements.begin(), elements.end()-1);
    return new Cigar(newElements);
}
