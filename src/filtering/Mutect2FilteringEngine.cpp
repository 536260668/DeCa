//
// Created by cluster on 22-11-5.
//

#include "Mutect2FilteringEngine.h"
#include "TumorEvidenceFilter.h"
#include "StrandArtifactFilter.h"

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

SomaticClusteringModel* Mutect2FilteringEngine::getSomaticClusteringModel() {
    return &somaticClusteringModel;
}

Mutect2FilteringEngine::Mutect2FilteringEngine(M2FiltersArgumentCollection &MTFAC, const std::string& normal) : somaticClusteringModel(MTFAC), normalSample(normal),
                                                                                                                thresholdCalculator(MTFAC.initialPosteriorThreshold, MTFAC.maxFalsePositiveRate, MTFAC.fScoreBeta){
    Mutect2VariantFilter * tmp = new TumorEvidenceFilter();
    filters.emplace_back(tmp);
    tmp = new StrandArtifactFilter();
    filters.emplace_back(tmp);
}

std::vector<int>
Mutect2FilteringEngine::sumStrandCountsOverSamples(const std::shared_ptr<VariantContext> &vc, bool includeTumor,
                                                   bool includeNormal) {
    std::vector<int> result = std::vector<int>(4, 0);
    std::shared_ptr<GenoTypesContext> genotype = vc->getGenotypes();
    std::vector<std::shared_ptr<Genotype>> * genotypes = genotype->getGenotypes();
    for(auto &gt : *genotypes) {
        if((includeTumor) && isTumor(gt.get()) || (includeNormal && isNormal(gt.get()))) {
            if(gt->hasExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY)) {
                std::vector<int> tmp = gt->getExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY).getAttributeAsIntVector();
                for(int i = 0; i < 4; i++) {
                    result[i] += tmp[i];
                }
            }
        }
    }
    return result;
}

Mutect2FilteringEngine::~Mutect2FilteringEngine() {
    for(auto filter : filters) {
        delete filter;
    }
}

void Mutect2FilteringEngine::accumulateData(const std::shared_ptr<VariantContext> &vc,
                                            std::shared_ptr<ReferenceContext> referenceContext) {
    bool flag = false;
    for(auto & allele : vc->getAlleles()) {
        if(allele->getIsNonReference() && !allele->getIsNonRefAllele()) {
            flag = true;
        }
    }
    if(!flag) {
        return;
    }
    ErrorProbabilities errorProbabilities = ErrorProbabilities(filters, vc, this, referenceContext);
    for(auto filter : filters) {
        filter->accumulateDataForLearning(vc, errorProbabilities, this);
    }
    std::vector<int> tumorADs = sumADsOverSamples(vc, true, false);
    std::vector<double> tumorLogOdds = Mutect2FilteringEngine::getTumorLogOdds(vc);
    somaticClusteringModel.record(tumorADs, tumorLogOdds, errorProbabilities.getTechnicalArtifactProbability(), errorProbabilities.getNonSomaticProbability(), vc);
    thresholdCalculator.addArtifactProbability(errorProbabilities.getErrorProbability());
}

