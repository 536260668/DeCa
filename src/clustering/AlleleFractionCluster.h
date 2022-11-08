//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_ALLELEFRACTIONCLUSTER_H
#define MUTECT2CPP_MASTER_ALLELEFRACTIONCLUSTER_H


#include "Datum.h"
#include <vector>
#include <string>

class AlleleFractionCluster {
    virtual double logLikelihood(Datum datum) = 0;

    virtual double logLikelihood(int totalCount, int altCount) = 0;

    virtual void learn(std::vector<Datum> data) = 0;

    virtual std::string toString() = 0;
};


#endif //MUTECT2CPP_MASTER_ALLELEFRACTIONCLUSTER_H
