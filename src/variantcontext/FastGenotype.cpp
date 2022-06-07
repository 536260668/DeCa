//
// Created by 梦想家xixi on 2021/11/30.
//

#include "FastGenotype.h"
#include <utility>

FastGenotype::FastGenotype(std::string sampleName, std::vector<std::shared_ptr<Allele>> &alleles, bool isPhased, int GQ, int DP,
                           int *AD, int ADLength, int *PL, int PLLength, const std::string& filters,
                           std::map<std::string, AttributeValue> extendedAttributes) : Genotype(std::move(sampleName), filters) , alleles(alleles), _isPhased(isPhased), GQ(GQ),
                           DP(DP), AD(AD), ADLength(ADLength), PL(PL), PLLength(PLLength), extendedAttributes(std::move(extendedAttributes)){}

int *FastGenotype::getAD(int &length) {
    length = ADLength;
    return AD;
}

int *FastGenotype::getPL(int &length) {
    length = PLLength;
    return PL;
}

FastGenotype::~FastGenotype() = default;
