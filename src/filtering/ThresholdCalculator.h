//
// Created by cluster on 22-11-10.
//

#ifndef MUTECT2CPP_MASTER_THRESHOLDCALCULATOR_H
#define MUTECT2CPP_MASTER_THRESHOLDCALCULATOR_H

#include <vector>

class ThresholdCalculator {
private:
    std::vector<double> artifactProbabilities;
    double maxFalseDiscoveryRate;
    double fScoreBeta;
    double threshold;

public:
    void addArtifactProbability(double artifactProbability);
    ThresholdCalculator(double initialThreshold, double maxFalseDiscoveryRate, double fScoreBeta);
};


#endif //MUTECT2CPP_MASTER_THRESHOLDCALCULATOR_H
