//
// Created by 梦想家xixi on 2021/11/1.
//

#ifndef MUTECT2CPP_MASTER_MATHUTILS_H
#define MUTECT2CPP_MASTER_MATHUTILS_H

#include "cache/DigammaCache.h"
#include "cache/Log10FactorialCache.h"


class MathUtils {
private:
    static DigammaCache DIGAMMA_CACHE;
    static Log10FactorialCache LOG_10_FACTORIAL_CACHE;
public:
    static double digamma(int n);

    static double log10ToLog(double log10);

    static double log10Factorial(int n);

    static double fastBernoulliEntropy(double p);
};


#endif //MUTECT2CPP_MASTER_MATHUTILS_H
