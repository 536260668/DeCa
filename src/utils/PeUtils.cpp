//
// Created by 梦想家xixi on 2022/1/4.
//

#include "PeUtils.h"
#include "Mutect2Utils.h"

PeUtils::PeUtils(SAMRecord *pe, int pos) : pos(pos), pe(pe), nCigarElements(&pe->getCigarElements()){
    bool flag = false;
    offset = 0;
    int start = pe->getStart();
    for(int i = 0; i < nCigarElements->size(); i++) {
        int length = (*nCigarElements)[i].getLength();
        CigarOperator tmp_cigarOperator = (*nCigarElements)[i].getOperator();
        if(!flag) {
            if(CigarOperatorUtils::getConsumesReferenceBases(tmp_cigarOperator)) {
                start += length;
            }
            if(CigarOperatorUtils::getConsumesReadBases(tmp_cigarOperator)) {
                offset += length;
            }
            if(start > pos){
                currentCigarElement = (*nCigarElements)[i];
                Cigar_offset = i;
                currentStart = (int)start - length;
                offset -=length;
                if(length != 0) {
                    offset = offset + pos - currentStart;
                }
                flag = true;
            }
        }

    }
    if(! flag)
        throw std::invalid_argument("pos and read is not match");
}

bool PeUtils::isBeforeSoftClip() {
    return isImmediatelyBefore(S);
}

bool PeUtils::isImmediatelyBefore(CigarOperator cigarOperator) {
    if(pos == currentStart + currentCigarElement.getLength() - 1) {
        if(Cigar_offset < (*nCigarElements).size() - 1 && cigarOperator == (*nCigarElements)[Cigar_offset+1].getOperator())
            return true;
    }
    return false;
}

bool PeUtils::isDeletion() {
    return currentCigarElement.getOperator() == D;
}

CigarElement &PeUtils::getCurrentCigarElement() {
    return currentCigarElement;
}

CigarElement* PeUtils::getNearestOnGenomeCigarElement(int direction) {
    for(int i = Cigar_offset + direction; i >= 0 && i < (*nCigarElements).size(); i += direction) {
        CigarElement* elt = &(*nCigarElements)[i];
        if(isOnGenomeCigar(elt->getOperator()))
            return elt;
    }
    return nullptr;
}

bool PeUtils::isOnGenomeCigar(CigarOperator cigarOperator) {
    if(cigarOperator == M || cigarOperator == EQ || cigarOperator == X || cigarOperator == D ) {
        return true;
    }
    return false;
}

int PeUtils::getLengthOfImmediatelyFollowingIndel() {
    CigarElement* element = getNextIndelCigarElement();
    return element != nullptr ? element->getLength() : 0;
}

CigarElement* PeUtils::getNextIndelCigarElement() {
    if(isBeforeSoftClip()) {
        return getNearestOnGenomeCigarElement(1);
    } else if(isBeforeInsertion()) {
        return &(*nCigarElements)[Cigar_offset + 1];
    } else
        return nullptr;
}

bool PeUtils::isBeforeDeletionStart() {
    if(currentCigarElement.getOperator() != D && currentStart + currentCigarElement.getLength() - 1 == pos
    && Cigar_offset < (*nCigarElements).size()-1 && (*nCigarElements)[Cigar_offset+1].getOperator() == D)
        return true;
    else
        return false;
}

bool PeUtils::isBeforeInsertion() {
    return isImmediatelyBefore(I);
}

uint8_t PeUtils::getQual() {
    if(isDeletion()) {
        return DELETION_QUAL;
    }
    return pe->getBaseQuality(offset);
}

uint8_t PeUtils::getBaseQuality(int pos) {
    if(isDeletion())
        return 'N';
    return pe->getBaseQuality(pos - pe->getStart());
}

bool PeUtils::isImmediatelyAfter(CigarOperator op) {
    return currentStart == pos && Cigar_offset - 1 >= 0 && (*nCigarElements)[Cigar_offset - 1].getOperator() == op;
}

bool PeUtils::isAfterSoftClip() {
    return isImmediatelyAfter(S);
}

uint8_t PeUtils::getBase() {
    Mutect2Utils::validateArg(pos <= pe->getEnd(), "pos is not in read");
    return pe->getBase(offset);
}
