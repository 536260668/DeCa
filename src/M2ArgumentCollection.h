//
// Created by lhh on 10/23/21.
//

#ifndef MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
#include <set>

typedef struct M2ArgumentCollection{
    int callableDepth;
    int maxProbPropagationDistance;
    double activeProbThreshold;
    int assemblyRegionPadding;
    int minAssemblyRegionSize;
    int maxAssemblyRegionSize;
    std::set<std::string> normalSamples;

}M2ArgumentCollection;

#endif //MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
