//
// Created by 梦想家xixi on 2021/11/9.
//

#ifndef MUTECT2CPP_MASTER_CIGAR_H
#define MUTECT2CPP_MASTER_CIGAR_H

#include "CigarElement.h"
#include <vector>

class Cigar {
private:
    std::vector<CigarElement> cigarElements;
    static bool isRealOperator(CigarOperator op);
    static bool isInDelOperator(CigarOperator op);
    static bool isClippingOperator(CigarOperator op);
    static bool isPaddingOperator(CigarOperator op);

public:
    Cigar()= default;
    explicit Cigar(std::vector<CigarElement>& cigarElements);
    std::vector<CigarElement> & getCigarElements();
    CigarElement & getCigarElement(int i);
    void add(CigarElement cigarElement);
    int numCigarElements() const {return (int)cigarElements.size();}
    bool isEmpty() const {return cigarElements.empty();}
    int getReferenceLength();
    int getPaddedReferenceLength();
    int getReadLength();
    static int getReadLength(std::vector<CigarElement> & cigarElements);
    static Cigar* fromCigarOperators(std::vector<CigarElement> & cigarElements);
    bool isLeftClipped();
    bool isRightClipped();
    bool isClipped();
    std::vector<std::logic_error> isValid(std::string readName, long recordNumber);
    CigarElement getFirstCigarElement();
    CigarElement getLastCigarElement();
};


#endif //MUTECT2CPP_MASTER_CIGAR_H
