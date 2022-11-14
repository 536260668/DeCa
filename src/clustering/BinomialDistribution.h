//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_BINOMIALDISTRIBUTION_H
#define MUTECT2CPP_MASTER_BINOMIALDISTRIBUTION_H


class BinomialDistribution {
private:
    int numberOfTrials;
    double probabilityOfSuccess;
    constexpr static double DEFAULT_EPSILON = 1E-14;
    static double getB(int n, double x, double a, double b);
    static double getA(int n, double x) {
        return 1.0;
    }

public:
    BinomialDistribution(int trails, double p);
    double cumulativeProbability(int x);
    static double regularizedBeta(double x, double a, double b);
    static double regularizedBeta(double x,
                                  double a, double b,
                                  double epsilon, int maxIterations);
    static double evaluate(double x, double epsilon, int maxIterations, double a, double b);
};


#endif //MUTECT2CPP_MASTER_BINOMIALDISTRIBUTION_H
