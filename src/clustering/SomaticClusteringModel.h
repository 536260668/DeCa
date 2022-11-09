//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H
#define MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H

#include "variantcontext/VariantContext.h"
#include "Datum.h"
#include "AlleleFractionCluster.h"
#include "BetaBinomialCluster.h"
#include "SequencingError.h"
#include "Mutect2FilteringEngine/M2FiltersArgumentCollection.h"


class SomaticClusteringModel {
private:
    std::map<int, double> logVariantPriors;
    std::vector<AlleleFractionCluster *> clusters;
    static BetaBinomialCluster NEW_CLUSTER;
    static const int SEQUENCING_ERROR_INDEX = 0;
    static const int HIGH_AF_INDEX = 1;
    static const int BACKGROUND_INDEX = 2;
    static const int OFFSET = 3;
    static const int MAX_INDEL_SIZE_IN_PRIOR_MAP = 10;
    constexpr static const double INITIAL_HIGH_AF_WEIGHT = 0.01;
    constexpr static const double INITIAL_BACKGROUND_WEIGHT = 0.01;
    constexpr static const double CONCENTRATION = 0.5;
    double logHighAFWeight;
    double logBackgroundWeight;
    double logSparseClustersWeight;
    static const BetaDistributionShape INITIAL_HIGH_AF_BETA;
    static const BetaDistributionShape INITIAL_BACKGROUND_BETA;

    std::vector<double > clusterProbabilities(Datum datum);
    std::vector<int> clusterCounts;
    double getLogPriorOfSomaticVariant(int indelLength);
    double logCRPWeight(int clusterIndex);
    int totalSparseClusterCount;
public:
    SomaticClusteringModel(M2FiltersArgumentCollection &MTFAC);
    double probabilityOfSequencingError(const Datum& datum);
    static int indelLength(const std::shared_ptr<VariantContext>& vc, int altIndex);
    virtual ~SomaticClusteringModel();
};


#endif //MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H
