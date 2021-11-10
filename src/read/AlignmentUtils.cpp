//
// Created by 梦想家xixi on 2021/11/10.
//

#include "AlignmentUtils.h"

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
        lastOp = cur.getOperator()
    }
    return false;
}
