//
// Created by cluster on 22-11-15.
//

#include "PolymeraseSlippageFilter.h"

PolymeraseSlippageFilter::PolymeraseSlippageFilter(int minSlippageLength, double slippageRate) : minSlippageLength(minSlippageLength), slippageRate(slippageRate){

}

ErrorType PolymeraseSlippageFilter::errorType() {
    return ARTIFACT;
}

double PolymeraseSlippageFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                           Mutect2FilteringEngine *filteringEngine,
                                                           std::shared_ptr<ReferenceContext>) {
    return 0;
}
