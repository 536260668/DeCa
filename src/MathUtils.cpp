//
// Created by lhh on 10/25/21.
//

#include <cmath>
#include "assert.h"
#include "MathUtils.h"
using namespace std;

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