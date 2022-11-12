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
#include "filtering/M2FiltersArgumentCollection.h"
#include "boost/random/mersenne_twister.hpp"
#include <random>

class SomaticClusteringModel {
private:
    std::map<int, double> logVariantPriors;
    std::vector<AlleleFractionCluster *> clusters;
    std::uniform_real_distribution<float> rng;
    static thread_local boost::random::mt19937 seed;
    static BetaBinomialCluster NEW_CLUSTER;
    static const int SEQUENCING_ERROR_INDEX = 0;
    static const int HIGH_AF_INDEX = 1;
    static const int BACKGROUND_INDEX = 2;
    static const int OFFSET = 3;
    static const int MAX_INDEL_SIZE_IN_PRIOR_MAP = 10;
    static const int NUM_ITERATIONS = 5;
    constexpr static const double INITIAL_HIGH_AF_WEIGHT = 0.01;
    constexpr static const double INITIAL_BACKGROUND_WEIGHT = 0.01;
    constexpr static const double CONCENTRATION = 0.5;
    double logHighAFWeight;
    double logBackgroundWeight;
    double logSparseClustersWeight;
    static BetaDistributionShape INITIAL_HIGH_AF_BETA;
    static BetaDistributionShape INITIAL_BACKGROUND_BETA;
    bool firstPass;

    std::vector<double > clusterProbabilities(Datum datum);
    std::vector<int> clusterCounts;
    std::vector<Datum> data;
    std::vector<std::optional<int>> clusterAssignments;
    double getLogPriorOfSomaticVariant(int indelLength);
    double logCRPWeight(int clusterIndex);
    SomaticClusteringModel(const SomaticClusteringModel& engine) = default;
    SomaticClusteringModel &operator=(const SomaticClusteringModel& engine) = default;
    Datum popDatum(int datumIndex);
    void pruneEmptyClusters();
    void learnWeightsAndPriors();
    double REGULARIZING_PSEUDOCOUNT;
    void assignDatum(int datumIndex, int clusterIndex);
    int totalSparseClusterCount;
public:
    SomaticClusteringModel(M2FiltersArgumentCollection &MTFAC);
    double probabilityOfSequencingError(const Datum& datum);
    static int indelLength(const std::shared_ptr<VariantContext>& vc, int altIndex);
    virtual ~SomaticClusteringModel();
    void record(const std::vector<int> & tumorADs, const std::vector<double> & tumorLogOdds, double artifactProbability, double nonSomaticProbability, const std::shared_ptr<VariantContext> & vc);
    void learnAndClearAccumulatedData();
};


#endif //MUTECT2CPP_MASTER_SOMATICCLUSTERINGMODEL_H
