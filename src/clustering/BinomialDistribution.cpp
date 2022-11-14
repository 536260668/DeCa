//
// Created by cluster on 22-11-14.
//

#include <climits>
#include <cmath>
#include "BinomialDistribution.h"
#include <stdexcept>

BinomialDistribution::BinomialDistribution(int trails, double p) : probabilityOfSuccess(p), numberOfTrials(trails){

}

double BinomialDistribution::cumulativeProbability(int x) {
    double ret;
    if (x < 0) {
        ret = 0.0;
    } else if (x >= numberOfTrials) {
        ret = 1.0;
    } else {
        ret = 1.0 - regularizedBeta(probabilityOfSuccess,
                                         x + 1.0, numberOfTrials - x);
    }
    return ret;
}

double BinomialDistribution::regularizedBeta(double x, double a, double b) {
    return regularizedBeta(x, a, b, DEFAULT_EPSILON, INT_MAX);
}

double BinomialDistribution::regularizedBeta(double x, double a, double b, double epsilon, int maxIterations) {
    double ret;

    if (isinf(x) ||
            isinf(a) ||
            isinf(b) ||
        x < 0 ||
        x > 1 ||
        a <= 0 ||
        b <= 0) {
        ret = nan("");
    } else if (x > (a + 1) / (2 + b + a) &&
               1 - x <= (b + 1) / (2 + b + a)) {
        ret = 1 - regularizedBeta(1 - x, b, a, epsilon, maxIterations);
    } else {

    }
}


double BinomialDistribution::evaluate(double x, double epsilon, int maxIterations, double a, double b) {
    double small = 1e-50;
    double hPrev = getA(0, x);

    // use the value of small as epsilon criteria for zero checks
    if (std::abs(hPrev) < small) {
        hPrev = small;
    }

    int n = 1;
    double dPrev = 0.0;
    double cPrev = hPrev;
    double hN = hPrev;

    while (n < maxIterations) {
        double a = getA(n, x);
        double b = getB(n, x, a, b);

        double dN = a + b * dPrev;
        if (std::abs(dN) < small) {
            dN = small;
        }
        double cN = a + b / cPrev;
        if (std::abs(cN) < small) {
            cN = small;
        }

        dN = 1 / dN;
        double deltaN = cN * dN;
        hN = hPrev * deltaN;

        if (std::isinf(hN)) {
            throw std::invalid_argument("");
        }
        if (std::isinf(hN)) {
            throw std::invalid_argument("");
        }

        if (std::abs(deltaN - 1.0) < epsilon) {
            break;
        }

        dPrev = dN;
        cPrev = cN;
        hPrev = hN;
        n++;
    }

    if (n >= maxIterations) {
        throw std::invalid_argument("");
    }

    return hN;
}

double BinomialDistribution::getB(int n, double x, double a, double b) {
    {
        double ret;
        double m;
        if (n % 2 == 0) { // even
            m = n / 2.0;
            ret = (m * (b - m) * x) /
                  ((a + (2 * m) - 1) * (a + (2 * m)));
        } else {
            m = (n - 1.0) / 2.0;
            ret = -((a + m) * (a + b + m) * x) /
                  ((a + (2 * m)) * (a + (2 * m) + 1.0));
        }
        return ret;
    }
}

