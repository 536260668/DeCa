//
// Created by cluster on 22-11-8.
//

#include "M2FiltersArgumentCollection.h"
#include "MathUtils.h"

const double M2FiltersArgumentCollection:: DEFAULT_LOG_INDEL_PRIOR = MathUtils::log10ToLog(-7);
const double M2FiltersArgumentCollection:: DEFAULT_LOG_SNV_PRIOR = MathUtils::log10ToLog(-6);
const double M2FiltersArgumentCollection:: DEFAULT_LOG_INDEL_PRIOR_FOR_MITO = MathUtils::log10ToLog(-3.75);
const double M2FiltersArgumentCollection:: DEFAULT_LOG_SNV_PRIOR_FOR_MITO = MathUtils::log10ToLog(-2.5);

double M2FiltersArgumentCollection::getLogIndelPrior() {
    return mitochondria && logIndelPrior == DEFAULT_LOG_INDEL_PRIOR ? DEFAULT_LOG_INDEL_PRIOR_FOR_MITO : logIndelPrior;
}

M2FiltersArgumentCollection::M2FiltersArgumentCollection() {
    mitochondria = false;
    logIndelPrior = DEFAULT_LOG_INDEL_PRIOR;
    logSNVPrior = DEFAULT_LOG_SNV_PRIOR;
}

double M2FiltersArgumentCollection::getLogSnvPrior() {
    return mitochondria && logSNVPrior == DEFAULT_LOG_SNV_PRIOR ? DEFAULT_LOG_SNV_PRIOR_FOR_MITO : logSNVPrior;
}
