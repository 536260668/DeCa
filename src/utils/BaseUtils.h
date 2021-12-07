//
// Created by 梦想家xixi on 2021/12/6.
//

#ifndef MUTECT2CPP_MASTER_BASEUTILS_H
#define MUTECT2CPP_MASTER_BASEUTILS_H


#include <cstdint>

class BaseUtils {
private:
    static int baseIndexMap[256];

public:
    void initial();
    static bool isRegularBase(uint8_t base);
    static int simpleBaseToBaseIndex(uint8_t base);
    static bool isAllRegularBases(uint8_t *bases, int length);
};


#endif //MUTECT2CPP_MASTER_BASEUTILS_H
