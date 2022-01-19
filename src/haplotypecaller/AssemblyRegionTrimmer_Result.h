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
    std::shared_ptr<AssemblyRegion> originalRegion;
    SimpleInterval* callableSpan;
    SimpleInterval* maximumSpan;
    SimpleInterval* extendedSpan;
    SimpleInterval* idealSpan;
    std::pair<SimpleInterval*, SimpleInterval*> * nonVariantFlanks;
    std::vector<std::shared_ptr<VariantContext>> callableEvents;
    int padding;
    int usableExtension;
    std::shared_ptr<AssemblyRegion> callableRegion;

private:
    std::shared_ptr<AssemblyRegion> leftFlankRegion;
    std::shared_ptr<AssemblyRegion> rightFlankRegion;
    bool emitReferenceConfidence;

public:
    AssemblyRegionTrimmer_Result(bool emitReferenceConfidence, bool needsTrimming, std::shared_ptr<AssemblyRegion> originalRegion, int padding, int extension,
                                 std::vector<std::shared_ptr<VariantContext>> * overlappingEvents, std::pair<SimpleInterval*, SimpleInterval*> * nonVariantFlanks,
                                 SimpleInterval* extendedSpan, SimpleInterval* idealSpan, SimpleInterval* maximumSpan, SimpleInterval* callableSpan);
    static std::shared_ptr<AssemblyRegionTrimmer_Result> noVariation(bool emitReferenceConfidence, std::shared_ptr<AssemblyRegion> targetRegion, int padding, int usableExtension);
    ~AssemblyRegionTrimmer_Result();
    static std::shared_ptr<AssemblyRegionTrimmer_Result> noTrimming(bool emitReferenceConfidence, std::shared_ptr<AssemblyRegion> targetRegion, int padding, int usableExtension, std::vector<std::shared_ptr<VariantContext>> * events);
    bool isVariationPresent();
    std::shared_ptr<AssemblyRegion> getCallableRegion();
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H
