//
// Created by cluster on 22-11-5.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2FILTERINGENGINE_H
#define MUTECT2CPP_MASTER_MUTECT2FILTERINGENGINE_H

#include "Mutect2VariantFilter.h"

class Mutect2FilteringEngine {
private:
    std::string normalSample;
public:
    bool isNormal(Genotype* genotype);
    bool isTumor(Genotype* genotype);
    static double roundFinitePrecisionErrors(double roundFinitePrecisionErrors);
    static std::vector<double> getTumorLogOdds(const std::shared_ptr<VariantContext> & vc);
    std::vector<int> sumADsOverSamples(const std::shared_ptr<VariantContext> & vc, bool includeTumor, bool includeNormal);
};


#endif //MUTECT2CPP_MASTER_MUTECT2FILTERINGENGINE_H
