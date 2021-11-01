/**
 * The implementation of QualityUtils class
 */
#include <iostream>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include "QualityUtils.h"

double QualityUtils::qualToErrorProb(int qual)
{
    return pow(10.0, ((double )qual) / -10.0);
}

char QualityUtils::errorProbToQual(double errorRate)
{
    return errorProbToQual(errorRate, MAX_SAM_QUAL_SCORE);
}

char QualityUtils::errorProbToQual(double errorRate, char maxQual)
{
    assert(errorRate >= 0.0 && errorRate <= 1.0);
    double d = round(-10.0 * log10(errorRate));
    int qual = isinf(d) ? INT32_MAX : (int)d;   // if d is infinity, (int)d will be -2147483647

    return boundQual(qual, maxQual);
}

char QualityUtils::boundQual(int qual, char maxQual)
{
    return (char)(std::max(std::min(qual, maxQual & 0xFF) ,1) & 0xFF);
}

