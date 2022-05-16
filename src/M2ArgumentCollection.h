//
// Created by lhh on 10/23/21.
//

#ifndef MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H

#include <set>
#include "haplotypecaller/LikelihoodEngineArgumentCollection.h"

typedef struct M2ArgumentCollection{
    int callableDepth;
    int maxProbPropagationDistance;
    double activeProbThreshold;
    int assemblyRegionPadding;
    int minAssemblyRegionSize;
    int maxAssemblyRegionSize;
    //std::set<std::string> normalSamples;
    std::string normalSample;
    bool genotypeGermlineSites = false;
    LikelihoodEngineArgumentCollection likelihoodArgs;
    int pcrSnvQual = 40;
    int pcrIndelQual = 40;

    static double getInitialLogOdds() {
		// A constant should not be returned, it caused a bug sometimes, one over a billion
		// So I won't implement MTAC class for the time being
		// 4.605170185988092 == 2 * ln(10)
		return 4.605170185988092;
    }
}M2ArgumentCollection;

#endif //MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
