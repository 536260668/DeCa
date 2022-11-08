//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_SEQUENCINGERROR_H
#define MUTECT2CPP_MASTER_SEQUENCINGERROR_H

#include "AlleleFractionCluster.h"

class SequencingError : AlleleFractionCluster{
public:
    double logLikelihood(Datum datum) override;

    double logLikelihood(int totalCount, int altCount) override;

    void learn(std::vector<Datum> data) override;

    std::string toString() override;
};


#endif //MUTECT2CPP_MASTER_SEQUENCINGERROR_H
