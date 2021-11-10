//
// Created by 梦想家xixi on 2021/11/1.
//

#include "MathUtils.h"
#include <cmath>

DigammaCache MathUtils::DIGAMMA_CACHE = DigammaCache();
Log10FactorialCache MathUtils::LOG_10_FACTORIAL_CACHE = Log10FactorialCache();

double MathUtils::digamma(int n) {
    return DIGAMMA_CACHE.get(n);
}

double MathUtils::log10ToLog(double log10) {
    return log10 * std::log(10);
}

double MathUtils::log10Factorial(int n) {
    return LOG_10_FACTORIAL_CACHE.get(n);
}

double MathUtils::fastBernoulliEntropy(const double p) {
    double product = p * (1 - p);
    return product * (11 + 33 * product) / (2 + 20 * product);
}
