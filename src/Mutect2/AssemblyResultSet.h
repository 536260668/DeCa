//
// Created by 梦想家xixi on 2021/12/1.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
#define MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H

#include <map>
#include "AssemblyResult.h"
#include "Haplotype.h"
#include "AssemblyRegion.h"
#include "VariantContext.h"

class AssemblyResultSet {
private:
    std::map<int, AssemblyResult*> assemblyResultByKmerSize;
    std::set<Haplotype*> haplotypes;
    std::map<Haplotype*, AssemblyResult*> assemblyResultByHaplotype;
    AssemblyRegion* regionForGenotyping;
    uint8_t * fullReferenceWithPadding;
    int fullReferenceWithPaddingLength;
    SimpleInterval* paddedReferenceLoc;
    bool variationPresent;
    Haplotype* refHaplotype;
    bool wasTrimmed = false;
    std::set<int> kmerSizes;
    std::set<VariantContext> variationEvents;
    bool add(AssemblyResult* ar);
    void updateReferenceHaplotype(Haplotype* newHaplotype);

public:
    AssemblyResultSet() = default;
    bool add(Haplotype * h, AssemblyResult* ar);
    bool add(Haplotype* h);
    void setRegionForGenotyping(AssemblyRegion* regionForGenotyping);
    void setFullReferenceWithPadding(uint8_t* fullReferenceWithPadding, int length);
    void setPaddedReferenceLoc(SimpleInterval* paddedReferenceLoc);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
