//
// Created by cluster on 22-11-9.
//

#ifndef MUTECT2CPP_MASTER_ERRORPROBABILITIES_H
#define MUTECT2CPP_MASTER_ERRORPROBABILITIES_H

#include "Mutect2VariantFilter.h"
#include <unordered_map>

class Mutect2VariantFilter;
class Mutect2FilteringEngine;

enum ErrorType {
    ARTIFACT, NON_SOMATIC, SEQUENCING
};


class ErrorProbabilities {
private:
    std::unordered_map<Mutect2VariantFilter *, double> probabilitiesByFilter;
    std::unordered_map<ErrorType, double> probabilitiesByType;
    double errorProbability;

public:
    ErrorProbabilities(std::vector<Mutect2VariantFilter *> filters, const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext);
    double getErrorProbability() const { return errorProbability; }
    double getTechnicalArtifactProbability() { return probabilitiesByType.at(ARTIFACT); }
    double getNonSomaticProbability() { return probabilitiesByType.at(NON_SOMATIC); }
    std::unordered_map<Mutect2VariantFilter *, double> getProbabilitiesByFilter() { return probabilitiesByFilter; }

};


#endif //MUTECT2CPP_MASTER_ERRORPROBABILITIES_H
