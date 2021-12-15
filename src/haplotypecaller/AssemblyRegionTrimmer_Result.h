//
// Created by 梦想家xixi on 2021/12/14.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H
#define MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H


#include "AssemblyRegion.h"
#include "VariantContext.h"

class AssemblyRegionTrimmer_Result {
protected:
    bool needsTrimming;
    AssemblyRegion* originalRegion;
    SimpleInterval* callableSpan;
    SimpleInterval* maximumSpan;
    SimpleInterval* extendedSpan;
    SimpleInterval* idealSpan;
    std::pair<SimpleInterval*, SimpleInterval*> * nonVariantFlanks;
    std::vector<VariantContext*> * callableEvents;
    int padding;
    int usableExtension;
    AssemblyRegion* callableRegion;

private:
    AssemblyRegion* leftFlankRegion;
    AssemblyRegion* rightFlankRegion;
    bool emitReferenceConfidence;

public:
    AssemblyRegionTrimmer_Result(bool emitReferenceConfidence, bool needsTrimming, AssemblyRegion* originalRegion, int padding, int extension,
                                 std::vector<VariantContext*> * overlappingEvents, std::pair<SimpleInterval*, SimpleInterval*> * nonVariantFlanks,
                                 SimpleInterval* extendedSpan, SimpleInterval* idealSpan, SimpleInterval* maximumSpan, SimpleInterval* callableSpan);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H
