//
// Created by 梦想家xixi on 2022/1/4.
//

#include "PeUtils.h"

PeUtils::PeUtils(bam1_t *pe, int pos) : pos(pos), pe(pe){
    uint32_t * res = bam_get_cigar(pe);
    hts_pos_t start = pe->core.pos;
    uint32_t n = pe->core.n_cigar;
    bool flag = false;
    for(int i = 0; i < n; i++) {
        int length = (int)(res[i] >> 4);
        CigarOperator tmp_cigarOperator = CigarOperatorUtils::binaryToEnum((int)(res[i] & 0xf));
        nCigarElements.emplace_back(CigarElement(length, tmp_cigarOperator));
        if(CigarOperatorUtils::getConsumesReferenceBases(tmp_cigarOperator)) {
            start += length;
        }
        if(start >= pos){
            currentCigarElement = CigarElement(length, tmp_cigarOperator);
            Cigar_offset = i;
            currentStart = (int)start - length;
            flag = true;
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
        if(Cigar_offset < nCigarElements.size() - 1 && cigarOperator == nCigarElements[Cigar_offset+1])
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

CigarElement PeUtils::getNearestOnGenomeCigarElement(int direction) {
    for(int i = Cigar_offset + direction; i >= 0 && i < nCigarElements.size(); i += direction) {
        CigarElement elt = nCigarElements[i];
        if(isOnGenomeCigar(elt.getOperator()))
            return elt;
    }
    return {0, M};
}

bool PeUtils::isOnGenomeCigar(CigarOperator cigarOperator) {
    if(cigarOperator == M || cigarOperator == EQ || cigarOperator == X || cigarOperator == D ) {
        return true;
    }
    return false;
}

int PeUtils::getLengthOfImmediatelyFollowingIndel() {
    CigarElement element = getNextIndelCigarElement();
    return element.getLength();
}

CigarElement PeUtils::getNextIndelCigarElement() {
    if(isBeforeSoftClip()) {
        return getNearestOnGenomeCigarElement(1);
    } else if(isBeforeInsertion()) {
        return nCigarElements[Cigar_offset + 1];
    } else
        return {0, M};
}

bool PeUtils::isBeforeDeletionStart() {
    if(currentCigarElement.getOperator() != D && currentStart + currentCigarElement.getLength() - 1 == pos
    && Cigar_offset < nCigarElements.size()-1 && nCigarElements[Cigar_offset+1].getOperator() == D)
        return true;
    else
        return false;
}

bool PeUtils::isBeforeInsertion() {
    return isImmediatelyBefore(I);
}
