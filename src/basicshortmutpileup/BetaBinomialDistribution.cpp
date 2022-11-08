//
// Created by cluster on 22-11-8.
//

#include "BetaBinomialDistribution.h"
#include "CombinatoricsUtils.h"
#include "boost/math/distributions.hpp"
#include <stdexcept>
#include <cfloat>

BetaBinomialDistribution::BetaBinomialDistribution(double alpha, double beta, int n, int rng) : alpha(alpha), beta(beta), n(n), rng(rng){
    if(alpha < 0 || beta < 0 || n <= 0) {
        std::string param = ("alpha, beta and n must be greater than zero.");
        throw std::invalid_argument(param);
    }
}

double BetaBinomialDistribution::logProbability(int k) {
    if(k < 0) {
        std::string param = ("alpha, beta and n must be greater than zero.");
        throw std::invalid_argument(param);
    }
    return k > n ? -DBL_MAX : CombinatoricsUtils::binomialCoefficientLog(n, k) + std::log(std::betaf(k+alpha, n - k + beta)) - std::log(std::betaf(alpha, beta));
}
