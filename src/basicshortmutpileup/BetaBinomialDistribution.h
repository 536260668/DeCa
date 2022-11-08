//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_BETABINOMIALDISTRIBUTION_H
#define MUTECT2CPP_MASTER_BETABINOMIALDISTRIBUTION_H


class BetaBinomialDistribution {
private:
    double alpha;
    double beta;
    int n;
    int rng;

public:
    BetaBinomialDistribution(double alpha, double beta, int n, int rng);
    double logProbability(int k);
};


#endif //MUTECT2CPP_MASTER_BETABINOMIALDISTRIBUTION_H
