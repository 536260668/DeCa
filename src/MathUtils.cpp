//
// Created by 梦想家xixi on 2021/11/1.
//

#include "assert.h"
#include "MathUtils.h"
#include <cmath>
using namespace std;



double MathUtils::digamma(int n) {
    return DIGAMMA_CACHE().get(n);
}

double MathUtils::log10ToLog(double log10) {
    return log10 * std::log(10);
}

double MathUtils::log10Factorial(int n) {
    return LOG_10_FACTORIAL_CACHE().get(n);
}

double MathUtils::fastBernoulliEntropy(const double p) {
    double product = p * (1 - p);
    return product * (11 + 33 * product) / (2 + 20 * product);
}

double MathUtils::normalDistribution(double mean, double sd, double x)
{
    assert(sd >= 0);
    assert((!isnan(mean)) && (!isinf(mean)) && (!isnan(sd)) && (!isinf(sd)) && !isnan(x) && (!isinf(x)));
    return exp(-(x - mean) * (x - mean) / (2.0 * sd * sd)) / (sd * sqrt(2.0 * M_PI));
}

std::vector<double> * MathUtils::normalizeSumToZero( std::vector<double> * array)
{
    assert(array != nullptr);
    if(array->size() == 0)
        return array;

    double sum = 0.0;
    for(auto &element: *array)
    {
        sum += element;
    }
    for(auto &element: *array)
    {
        element = element/sum;
    }
    return array;
}

DigammaCache &MathUtils::DIGAMMA_CACHE() {
    static DigammaCache digamma;
    return digamma;
}

Log10FactorialCache &MathUtils::LOG_10_FACTORIAL_CACHE() {
    static Log10FactorialCache log10Factorial;
    return log10Factorial;
}


