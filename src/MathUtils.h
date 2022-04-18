/**
 * A helper class to math calculation
 */

#ifndef MUTECT2CPP_MASTER_MATHUTILS_H
#define MUTECT2CPP_MASTER_MATHUTILS_H

#include <vector>
#include "cache/DigammaCache.h"
#include "cache/Log10FactorialCache.h"
#include "AssemblyRegion.h"
#include "ReadCache.h"


class MathUtils {
private:
    static DigammaCache & DIGAMMA_CACHE();
    static Log10FactorialCache & LOG_10_FACTORIAL_CACHE();
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

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    static int fastRound(double d);

    /*
     * bionomial Probability(int, int, double ) with log applied to result
     */
    static double log10BinomialProbability(int n, int k, double log10p);

    /*
     *  Calculates the log10 of the gamma function for x
     * @param x
     * @return
     */
    static double log10Gamma(double x);

    static double log10Fractorial(int n);

    static double log10BinomialCoefficient(int n, int k);
};


#endif //MUTECT2CPP_MASTER_MATHUTILS_H
