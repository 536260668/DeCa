//
// Created by cluster on 22-11-8.
//

#include "SomaticClusteringModel.h"
#include "NaturalLogUtils.h"

BetaBinomialCluster SomaticClusteringModel::NEW_CLUSTER(BetaDistributionShape::FLAT_BETA);
const BetaDistributionShape SomaticClusteringModel::INITIAL_HIGH_AF_BETA = BetaDistributionShape(10, 1);
const BetaDistributionShape SomaticClusteringModel::INITIAL_BACKGROUND_BETA = BetaDistributionShape::FLAT_BETA;

std::vector<double> SomaticClusteringModel::clusterProbabilities(Datum datum) {
    double logVariantPrior = getLogPriorOfSomaticVariant(datum.getIndelLength());
    double logNoVariantPrior = NaturalLogUtils::log1mexp(logVariantPrior);

    std::vector<double> logClusterPosteriors = std::vector<double>(clusters.size() + 1, 0);
    for(int i = 0; i < clusters.size() + 1; i++) {
        double logLikelihood = i < clusters.size() ? clusters[i]->logLikelihood(datum) : NEW_CLUSTER.logLikelihood(datum);
        if(i == SEQUENCING_ERROR_INDEX) {
            logClusterPosteriors[i] = logNoVariantPrior + logLikelihood;
        } else if (i == HIGH_AF_INDEX) {
            logClusterPosteriors[i] =  logVariantPrior + logHighAFWeight + logLikelihood;
        } else if (i == BACKGROUND_INDEX) {
            logClusterPosteriors[i] = logVariantPrior + logBackgroundWeight + logLikelihood;
        } else if (i < clusters.size()) {   // existing sparse cluster
            logClusterPosteriors[i] = logVariantPrior + logSparseClustersWeight + logCRPWeight(i)
                   + logLikelihood;
        } else {    // new sparse cluster
            logClusterPosteriors[i] = logVariantPrior + logSparseClustersWeight + logCRPWeight(i)
                   + logLikelihood;
        }
    }

    return  NaturalLogUtils::normalizeLog(logClusterPosteriors, false, false);
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

SomaticClusteringModel::SomaticClusteringModel(M2FiltersArgumentCollection &MTFAC) {
    totalSparseClusterCount = 0;
    logHighAFWeight = std::log(INITIAL_HIGH_AF_WEIGHT);
    logBackgroundWeight = std::log(INITIAL_BACKGROUND_WEIGHT);
    std::vector<double> inputs = {logHighAFWeight, logBackgroundWeight};
    logSparseClustersWeight = NaturalLogUtils::log1mexp(NaturalLogUtils::logSumExp(inputs));
    for(int i = -MAX_INDEL_SIZE_IN_PRIOR_MAP; i < MAX_INDEL_SIZE_IN_PRIOR_MAP+1; i++) {
        logVariantPriors.insert({i, MTFAC.getLogIndelPrior()});
    }
    logVariantPriors[0] = MTFAC.getLogSnvPrior();
    clusters = std::vector<AlleleFractionCluster *>(3);
    clusters[SEQUENCING_ERROR_INDEX] = new SequencingError();
    clusters[HIGH_AF_INDEX] = new BetaBinomialCluster(INITIAL_HIGH_AF_BETA);
    clusters[BACKGROUND_INDEX] = new BetaBinomialCluster(INITIAL_BACKGROUND_BETA);
}

double SomaticClusteringModel::logCRPWeight(int clusterIndex) {
    if(clusterIndex >= OFFSET) {
        throw std::invalid_argument("Chinese restaurant process does not apply to error, high-AF, and backgorund clusters");
    }
    double numerator = clusterIndex == clusters.size() ? CONCENTRATION : clusterCounts[clusterIndex];
    double denominator = totalSparseClusterCount + CONCENTRATION;
    return std::log(numerator/ denominator);
}

double SomaticClusteringModel::probabilityOfSequencingError(const Datum &datum) {
    return clusterProbabilities(datum)[SEQUENCING_ERROR_INDEX];
}

int SomaticClusteringModel::indelLength(const std::shared_ptr<VariantContext> &vc, int altIndex) {
    return vc->getAlternateAllele(altIndex)->getLength() - vc->getReference()->getLength();
}

SomaticClusteringModel::~SomaticClusteringModel() {
    for(auto c : clusters) {
        delete c;
    }
}
