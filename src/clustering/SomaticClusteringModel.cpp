//
// Created by cluster on 22-11-8.
//

#include "SomaticClusteringModel.h"
#include "NaturalLogUtils.h"

std::vector<double> SomaticClusteringModel::clusterProbabilities(Datum datum) {
    double logVariantPrior = getLogPriorOfSomaticVariant(datum.getIndelLength());
    double logNoVariantPrior = NaturalLogUtils::log1mexp(logVariantPrior);

}

double SomaticClusteringModel::getLogPriorOfSomaticVariant(int indelLength) {
    if(logVariantPriors.find(indelLength) == logVariantPriors.end()) {
        double input = (*logVariantPriors.begin()).second;
        for(auto & kv : logVariantPriors) {
            if(kv.second < input) {
                input = kv.second;
            }
        }
        logVariantPriors.insert({indelLength, input});
    }
    return logVariantPriors[indelLength] + (indelLength == 0 ? std::log(3.0) : 0);
}
