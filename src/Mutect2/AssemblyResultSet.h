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
//#include "ReadThreadingAssembler.h"

struct HaplotypeComp
{
public:
    bool operator()(const std::shared_ptr<Haplotype>& left, const std::shared_ptr<Haplotype>& right)
    {
        return (*left) < (*right);
    }
};

class AssemblyResultSet {
private:
    std::map<int, std::shared_ptr<AssemblyResult>> assemblyResultByKmerSize;
    std::set<std::shared_ptr<Haplotype>, HaplotypeComp> haplotypes;
    std::map<std::shared_ptr<Haplotype>, std::shared_ptr<AssemblyResult>> assemblyResultByHaplotype;
    std::shared_ptr<AssemblyRegion> regionForGenotyping;
    std::shared_ptr<uint8_t[]> fullReferenceWithPadding;
    int fullReferenceWithPaddingLength;
    SimpleInterval paddedReferenceLoc;
    bool variationPresent;
    std::shared_ptr<Haplotype>  refHaplotype;
    bool wasTrimmed = false;
    int lastMaxMnpDistanceUsed = -1;
    std::set<int> kmerSizes;
    std::set<std::shared_ptr<VariantContext>, VariantContextComparator> variationEvents;
    bool add(std::shared_ptr<AssemblyResult> &ar);
    void updateReferenceHaplotype(std::shared_ptr<Haplotype> & newHaplotype);

public:
    AssemblyResultSet() = default;
    bool add(std::shared_ptr<Haplotype> & h, std::shared_ptr<AssemblyResult> &ar);
    bool add(std::shared_ptr<Haplotype> & h);
    void setRegionForGenotyping(std::shared_ptr<AssemblyRegion> regionForGenotyping);
    void setFullReferenceWithPadding(std::shared_ptr<uint8_t[]> fullReferenceWithPadding, int length);
    void setPaddedReferenceLoc(SimpleInterval* paddedReferenceLoc);
    std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & getVariationEvents(int maxMnpDistance);

    void regenerateVariationEvents(int distance);
    std::vector<std::shared_ptr<Haplotype>> getHaplotypeList();
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
