//
// Created by cluster on 22-11-5.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2VARIANTFILTER_H
#define MUTECT2CPP_MASTER_MUTECT2VARIANTFILTER_H

#include "variantcontext/VariantContext.h"
#include "Mutect2FilteringEngine.h"
#include "engine/ReferenceContext.h"
#include "ErrorProbabilities.h"

class Mutect2FilteringEngine;
class ErrorProbabilities;

class Mutect2VariantFilter {
public:
    Mutect2VariantFilter() = default;
    double errorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine filteringEngine, const std::shared_ptr<ReferenceContext>& referenceContext);
    virtual ~Mutect2VariantFilter() {}

protected:
    virtual std::vector<std::string> requiredAnnotations() = 0;
    virtual double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine filteringEngine, std::shared_ptr<ReferenceContext>) = 0;
    virtual void learnParameters();
    virtual void learnParametersAndClearAccumulatedData();
    virtual void accumulateDataForLearning(const std::shared_ptr<VariantContext> & vc, ErrorProbabilities errorProbabilities, const Mutect2FilteringEngine& filteringEngine);
    virtual void clearAccumulatedData();
};


#endif //MUTECT2CPP_MASTER_MUTECT2VARIANTFILTER_H
