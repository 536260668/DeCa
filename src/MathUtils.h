//
// Created by 梦想家xixi on 2021/11/1.
//

#ifndef MUTECT2CPP_MASTER_MATHUTILS_H
#define MUTECT2CPP_MASTER_MATHUTILS_H

#include <vector>
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

    /**
     * Calculate f(x) = Normal(x | mu = mean, sigma = sd)
     * @param mean the desired mean of the Normal distribution
     * @param sd the desired standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    static double normalDistribution(double mean, double sd, double x);

    static std::vector<double> * normalizeSumToZero(std::vector<double> * array);
};


#endif //MUTECT2CPP_MASTER_MATHUTILS_H
