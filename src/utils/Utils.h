//
// Created by lhh on 4/26/22.
//

#ifndef MUTECT2CPP_MASTER_UTILS_H
#define MUTECT2CPP_MASTER_UTILS_H

#include <cstdint>
#include <memory>

class Utils {
public:
    static bool equalRange(uint8_t* left, int leftOffset, uint8_t* right, int rightOffset, int length);

    static std::shared_ptr<char[]> dupBytes(char b, int nCopies);
};


#endif //MUTECT2CPP_MASTER_UTILS_H
