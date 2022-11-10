//
// Created by cluster on 22-11-10.
//

#include "ThresholdCalculator.h"

void ThresholdCalculator::addArtifactProbability(double artifactProbability) {
    artifactProbabilities.emplace_back(artifactProbability);
}

ThresholdCalculator::ThresholdCalculator(double initialThreshold, double maxFalseDiscoveryRate, double fScoreBeta) : threshold(initialThreshold), maxFalseDiscoveryRate(maxFalseDiscoveryRate), fScoreBeta(fScoreBeta){

}
