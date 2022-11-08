//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_BETABINOMIALCLUSTER_H
#define MUTECT2CPP_MASTER_BETABINOMIALCLUSTER_H

#include "AlleleFractionCluster.h"
#include "BetaDistributionShape.h"


class BetaBinomialCluster : AlleleFractionCluster{
private:
    static const double RATE;
    static const int NUM_EPOCHS;
    BetaDistributionShape betaDistributionShape;

    static double logOddsCorrection(const BetaDistributionShape& originalBeta, const BetaDistributionShape& newBeta, int altCount, int refCount);

public:
    BetaBinomialCluster(const BetaDistributionShape& betaDistributionShape);

    static double logLikelihood(const Datum& datum, const BetaDistributionShape& betaDistributionShape);

    double logLikelihood(Datum datum) override;

    double logLikelihood(int totalCount, int altCount) override;

    void learn(std::vector<Datum> data) override;

    std::string toString() override;
};


#endif //MUTECT2CPP_MASTER_BETABINOMIALCLUSTER_H
