//
// Created by 梦想家xixi on 2021/10/30.
//

#ifndef MUTECT2CPP_MASTER_QUALITYUTILS_H
#define MUTECT2CPP_MASTER_QUALITYUTILS_H


#include <cstdint>

class QualityUtils {
private:
    static double qualToErrorProbCache[255];
    static void initial();

public:
    static uint8_t errorProbToQual(double prob, uint8_t maxQual);
    static uint8_t errorProbToQual(double errorRate);
    static uint8_t boundQual(int qual, uint8_t maxQual);
    static double qualToErrorProb(double qual);
    static double qualToErrorProb(uint8_t qual);
};


#endif //MUTECT2CPP_MASTER_QUALITYUTILS_H
