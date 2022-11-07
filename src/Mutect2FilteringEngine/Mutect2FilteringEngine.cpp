//
// Created by cluster on 22-11-5.
//

#include "Mutect2FilteringEngine.h"

double Mutect2FilteringEngine::roundFinitePrecisionErrors(double probability) {
    return std::max(std::min(probability, 1.0), 0.0);
}

std::vector<double> Mutect2FilteringEngine::getTumorLogOdds(const std::shared_ptr<VariantContext> &vc) {
    if(vc->hasAttribute(VCFConstants::TUMOR_LOG_10_ODDS_KEY)) {
        auto atrs = vc->getAttributes();
        AttributeValue atr = atrs.at(VCFConstants::TUMOR_LOG_10_ODDS_KEY);
        return atr.getAttributeAsDoubleVector();
    } else {
        return {};
    }
}

std::vector<int>
Mutect2FilteringEngine::sumADsOverSamples(const std::shared_ptr<VariantContext> &vc, bool includeTumor, bool includeNormal) {
    std::vector<int> ADs = std::vector<int>(vc->getNAlleles(), 0);
    std::shared_ptr<GenoTypesContext> genotype = vc->getGenotypes();
    std::vector<std::shared_ptr<Genotype>> * genotypes = genotype->getGenotypes();
    for(auto &gt : *genotypes) {
        if((includeTumor) && isTumor(gt.get()) || (includeNormal && isNormal(gt.get()))) {
            int len;
            int * ads = gt->getAD(len);
            for(int i = 0; i < vc->getNAlleles(); i++) {
                ADs[i] += ads[i];
            }
        }
    }
    return ADs;
}

bool Mutect2FilteringEngine::isNormal(Genotype *genotype) {
    std::string gt_name = genotype->getSampleName();
    if(normalSample.find(gt_name) != std::string::npos) {
        return true;
    }
    return false;
}

bool Mutect2FilteringEngine::isTumor(Genotype *genotype) {
    return !isNormal(genotype);
}
