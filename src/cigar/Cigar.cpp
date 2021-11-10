//
// Created by 梦想家xixi on 2021/11/9.
//

#include "Cigar.h"

Cigar::Cigar(std::vector<CigarElement>& cigarElements) {
    for(CigarElement e : cigarElements)
        this->cigarElements.emplace_back(e);
}

std::vector<CigarElement> & Cigar::getCigarElements() {
    return cigarElements;
}

CigarElement &Cigar::getCigarElement(int i) {
    return cigarElements[i];
}

void Cigar::add(CigarElement cigarElement) {
    cigarElements.emplace_back(cigarElement);
}

int Cigar::getReferenceLength() {
    int length = 0;
    for(CigarElement e : cigarElements) {
        switch (e.getOperator()) {
            case M:
            case D:
            case N:
            case EQ:
            case X:
                length += e.getLength();
        }
    }
    return length;
}

int Cigar::getPaddedReferenceLength() {
    int length = 0;
    for(CigarElement e : cigarElements) {
        switch (e.getOperator()) {
            case M:
            case D:
            case N:
            case EQ:
            case X:
            case P:
                length += e.getLength();
        }
    }
    return length;
}

int Cigar::getReadLength(std::vector<CigarElement> & cigarElements) {
    int length = 0;
    for(CigarElement e : cigarElements) {
        if(CigarOperatorUtils::getConsumesReadBases(e.getOperator())) {
            length += e.getLength();
        }
    }

    return length;
}

int Cigar::getReadLength() {
    return getReadLength(cigarElements);
}

Cigar* Cigar::fromCigarOperators(std::vector<CigarElement> &cigarElements) {
    if(cigarElements.empty())
        throw std::invalid_argument("cigarOperators is null");
    else {
        std::vector<CigarElement> cigarElementList;
        int j;
        for(int i = 0; i < cigarElements.size(); i = j) {
            CigarOperator currentOp = cigarElements[i].getOperator();
            for(j = i + 1; j < cigarElements.size() && cigarElements[j].getOperator() == currentOp; ++j){}
            cigarElementList.emplace_back(CigarElement(j - i, currentOp));
        }
        Cigar* res = new Cigar(cigarElementList);
        return res;
    }
}

bool Cigar::isRealOperator(CigarOperator op) {
    return op == M || op == EQ || op == X || op == I || op == D || op == N;
}

bool Cigar::isInDelOperator(CigarOperator op) {
    return CigarOperatorUtils::isIndel(op);
}

bool Cigar::isClippingOperator(CigarOperator op) {
    return CigarOperatorUtils::isClipping(op);
}

bool Cigar::isPaddingOperator(CigarOperator op) {
    return CigarOperatorUtils::isPadding(op);
}

bool Cigar::isLeftClipped() {
    return !isEmpty() && isClippingOperator(cigarElements[0].getOperator());
}

bool Cigar::isRightClipped() {
    return !isEmpty() && isClippingOperator(cigarElements[numCigarElements()-1].getOperator());
}

bool Cigar::isClipped() {
    return isLeftClipped() || isRightClipped();
}
