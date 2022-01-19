//
// Created by 梦想家xixi on 2022/1/15.
//

#ifndef MUTECT2CPP_MASTER_REFERENCECONFIDENCEMODEL_H
#define MUTECT2CPP_MASTER_REFERENCECONFIDENCEMODEL_H

#include "Haplotype.h"
#include "AssemblyRegion.h"

class ReferenceConfidenceModel {
public:
    static std::shared_ptr<Haplotype> createReferenceHaplotype(AssemblyRegion & activeRegion, uint8_t* refBase, int &length,SimpleInterval& paddedReferenceLoc);
};


#endif //MUTECT2CPP_MASTER_REFERENCECONFIDENCEMODEL_H
