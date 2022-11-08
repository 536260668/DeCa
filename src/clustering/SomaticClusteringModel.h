//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H
#define MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H

#include "variantcontext/VariantContext.h"
#include "Datum.h"
#include "AlleleFractionCluster.h"


class SomaticClusteringModel {
private:
    std::map<int, double> logVariantPriors;
    std::vector<AlleleFractionCluster> clusters;
    std::vector<double > clusterProbabilities(Datum datum);
    double getLogPriorOfSomaticVariant(int indelLength);

};


#endif //MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H
