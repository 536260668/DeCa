//
// Created by cluster on 22-11-5.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2FILTERINGENGINE_H
#define MUTECT2CPP_MASTER_MUTECT2FILTERINGENGINE_H

#include "Mutect2VariantFilter.h"
#include "SomaticClusteringModel.h"
#include "ErrorProbabilities.h"
#include "engine/ReferenceContext.h"
#include "ThresholdCalculator.h"

class SomaticClusteringModel;
class Mutect2VariantFilter;

class Mutect2FilteringEngine {
private:
    std::string normalSample;
    SomaticClusteringModel somaticClusteringModel;
    ThresholdCalculator thresholdCalculator;
    std::vector<Mutect2VariantFilter *> filters;
    Mutect2FilteringEngine(const Mutect2FilteringEngine& engine) = default;
    Mutect2FilteringEngine &operator=(const Mutect2FilteringEngine& engine) = default;

public:
    Mutect2FilteringEngine(M2FiltersArgumentCollection& MTFAC, const std::string& normal);
    bool isNormal(Genotype* genotype);
    bool isTumor(Genotype* genotype);
    static double roundFinitePrecisionErrors(double roundFinitePrecisionErrors);
    static std::vector<double> getTumorLogOdds(const std::shared_ptr<VariantContext> & vc);
    std::vector<int> sumADsOverSamples(const std::shared_ptr<VariantContext> & vc, bool includeTumor, bool includeNormal);
    SomaticClusteringModel* getSomaticClusteringModel();
    std::vector<int> sumStrandCountsOverSamples(const std::shared_ptr<VariantContext> & vc, bool includeTumor, bool includeNormal);
    ~Mutect2FilteringEngine();
    void accumulateData(const std::shared_ptr<VariantContext> & vc, std::shared_ptr<ReferenceContext> referenceContext);
    void learnParameters();
};


#endif //MUTECT2CPP_MASTER_MUTECT2FILTERINGENGINE_H
