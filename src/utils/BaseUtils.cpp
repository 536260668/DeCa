//
// Created by 梦想家xixi on 2021/12/6.
//

#include "BaseUtils.h"
#include "Mutect2Utils.h"

int BaseUtils::baseIndexMap[256] = {-1};

bool BaseUtils::isRegularBase(const uint8_t base) {
    return simpleBaseToBaseIndex(base) != -1;
}

int BaseUtils::simpleBaseToBaseIndex(uint8_t base) {
    Mutect2Utils::validateArg(base >= 0 && base < 256, "Non-standard bases were encountered in either the input reference or BAM file(s)");
    return baseIndexMap[base];
}

void BaseUtils::initial() {
    baseIndexMap['A'] = 0;
    baseIndexMap['a'] = 0;
    baseIndexMap['*'] = 0;    // the wildcard character counts as an A
    baseIndexMap['C'] = 1;
    baseIndexMap['c'] = 1;
    baseIndexMap['G'] = 2;
    baseIndexMap['g'] = 2;
    baseIndexMap['T'] = 3;
    baseIndexMap['t'] = 3;
}

bool BaseUtils::isAllRegularBases(std::shared_ptr<uint8_t[]> bases_, const int length) {
    uint8_t * bases = bases_.get();
    for(int i = 0; i < length; i++) {
        if(!isRegularBase(bases[i]))
            return false;
    }
    return true;
}
