//
// Created by cluster on 22-11-5.
//

#include "Mutect2VariantFilter.h"

double Mutect2VariantFilter::errorProbability(const std::shared_ptr<VariantContext> &vc,
                                              Mutect2FilteringEngine filteringEngine,
                                              const std::shared_ptr<ReferenceContext> &referenceContext) {
    bool flag = true;
    for(const auto & str : requiredAnnotations()) {
        if(!vc->hasAttribute(str)) {
            flag = false;
            break;
        }
    }
    double result = flag ? calculateErrorProbability(vc, filteringEngine, referenceContext): 0;
    return Mutect2FilteringEngine::roundFinitePrecisionErrors(result);
}
