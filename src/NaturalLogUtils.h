//
// Created by 梦想家xixi on 2021/11/1.
//

#ifndef MUTECT2CPP_MASTER_NATURALLOGUTILS_H
#define MUTECT2CPP_MASTER_NATURALLOGUTILS_H

#include <cmath>
#include <cstdint>

class NaturalLogUtils {
private:
    static constexpr double PHRED_TO_LOG_ERROR_PROB_FACTOR = -0.23025850929940458;
    static constexpr double LOG1MEXP_THRESHOLD = -0.6931471805599453;
    static double qualToLogProbCache[255];

public:
    // initialization of qualToLogProbCache[], you need to call this method before using any method of this class
    static void initial();
    static double qualToLogErrorProb(double qual);
    static double qualToLogErrorProb(uint8_t qual);
    static double log1mexp(double a);
    static double qualToLogProb(uint8_t qual);
};


#endif //MUTECT2CPP_MASTER_NATURALLOGUTILS_H
