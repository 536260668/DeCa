//
// Created by lhh on 10/25/21.
//

#ifndef MUTECT2CPP_MASTER_MATHUTILS_H
#define MUTECT2CPP_MASTER_MATHUTILS_H

#include <vector>

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful
 */
class MathUtils {
private:
    //constexpr static double ROOT_TWO_PI = sqrt(2.0 * M_PI);
public:
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
