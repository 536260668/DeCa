//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_DANGLINGCHAINMERGEHELPER_H
#define MUTECT2CPP_MASTER_DANGLINGCHAINMERGEHELPER_H

#include "MultiDeBruijnVertex.h"
#include "Cigar.h"
#include <vector>
#include <utility>

class DanglingChainMergeHelper {
public:
    std::vector<MultiDeBruijnVertex*> danglingPath;
    std::vector<MultiDeBruijnVertex*> referencePath;
    uint8_t * danglingPathString;
    uint8_t * referencePathString;
    int danglingPathStringLength;
    int referencePathStringLength;
    std::shared_ptr<Cigar> cigar;

    DanglingChainMergeHelper(std::vector<MultiDeBruijnVertex*> danglingPath, std::vector<MultiDeBruijnVertex*> referencePath, uint8_t * danglingPathString, int danglingPathStringLength,
                             uint8_t * referencePathString,  int referencePathStringLength, std::shared_ptr<Cigar> cigar) : danglingPath(std::move(danglingPath)), referencePath(std::move(referencePath)), danglingPathString(danglingPathString),
                             referencePathString(referencePathString), danglingPathStringLength(danglingPathStringLength), referencePathStringLength(referencePathStringLength), cigar(cigar) {}
};


#endif //MUTECT2CPP_MASTER_DANGLINGCHAINMERGEHELPER_H
