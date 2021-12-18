//
// Created by 梦想家xixi on 2021/12/14.
//

#include "AssemblyRegionTrimmer_Result.h"

AssemblyRegionTrimmer_Result::AssemblyRegionTrimmer_Result(bool emitReferenceConfidence, bool needsTrimming,
                                                           AssemblyRegion *originalRegion, int padding, int extension,
                                                           std::vector<VariantContext *> * overlappingEvents,
                                                           std::pair<SimpleInterval* , SimpleInterval* > * nonVariantFlanks,
                                                           SimpleInterval *extendedSpan, SimpleInterval *idealSpan,
                                                           SimpleInterval *maximumSpan, SimpleInterval *callableSpan) : emitReferenceConfidence(emitReferenceConfidence), needsTrimming(needsTrimming), callableEvents(*overlappingEvents),padding(padding),
                                                                                                                        usableExtension(extension)
                                                          {
    Mutect2Utils::validateArg(extendedSpan == nullptr || callableSpan == nullptr || extendedSpan->contains(callableSpan), "the extended callable span must include the callable span");
    if(originalRegion != nullptr)
        this->originalRegion = originalRegion;
    else
        this->originalRegion = nullptr;
    if(nonVariantFlanks != nullptr){
        this->nonVariantFlanks = new std::pair<SimpleInterval*, SimpleInterval*>();
        if(nonVariantFlanks->first != nullptr)
            this->nonVariantFlanks->first = new SimpleInterval(nonVariantFlanks->first);
        if(nonVariantFlanks->second != nullptr)
            this->nonVariantFlanks->second = new SimpleInterval(nonVariantFlanks->second);
    }
    else
        this->nonVariantFlanks = nullptr;
    if(extendedSpan != nullptr)
        this->extendedSpan = new SimpleInterval(extendedSpan);
    else
        this->extendedSpan = nullptr;
    if(idealSpan != nullptr)
        this->idealSpan = new SimpleInterval(idealSpan);
    else
        this->idealSpan = nullptr;
    if(maximumSpan != nullptr)
        this->maximumSpan = new SimpleInterval(maximumSpan);
    else
        this->maximumSpan = nullptr;
    if(callableSpan != nullptr)
        this->callableSpan = new SimpleInterval(callableSpan);
    else
        this->callableSpan = nullptr;

}

AssemblyRegionTrimmer_Result *
AssemblyRegionTrimmer_Result::noVariation(bool emitReferenceConfidence, AssemblyRegion *targetRegion, int padding,
                                          int usableExtension) {
    std::vector<VariantContext *> events;
    std::pair<SimpleInterval* , SimpleInterval* > nonVariantFlanks = std::pair<SimpleInterval*, SimpleInterval*>(&targetRegion->getSpan(), nullptr);
    AssemblyRegionTrimmer_Result* result = new AssemblyRegionTrimmer_Result(emitReferenceConfidence, false, targetRegion, padding, usableExtension,
                                                                            &events, &nonVariantFlanks,
                                                                            nullptr, nullptr, nullptr, nullptr);
    result->leftFlankRegion = targetRegion;
    return result;
}

AssemblyRegionTrimmer_Result *
AssemblyRegionTrimmer_Result::noTrimming(bool emitReferenceConfidence, AssemblyRegion *targetRegion, int padding,
                                         int usableExtension, std::vector<VariantContext *> *events) {
    SimpleInterval& targetRegionLoc = targetRegion->getSpan();
    std::pair<SimpleInterval* , SimpleInterval* > nonVariantFlanks = std::pair<SimpleInterval*, SimpleInterval*>(nullptr, nullptr);
    AssemblyRegionTrimmer_Result* result = new AssemblyRegionTrimmer_Result(emitReferenceConfidence, false, targetRegion, padding, usableExtension, events, &nonVariantFlanks, &targetRegionLoc, &targetRegionLoc, &targetRegionLoc, &targetRegionLoc);
    result->callableRegion = targetRegion;
    return result;
}

AssemblyRegionTrimmer_Result::~AssemblyRegionTrimmer_Result() {
    if(originalRegion == callableRegion)
        callableRegion = nullptr;
    if(originalRegion == leftFlankRegion)
        callableRegion = nullptr;
    if(originalRegion == rightFlankRegion)
        rightFlankRegion = nullptr;
    delete originalRegion;
    delete callableRegion;
    delete leftFlankRegion;
    delete rightFlankRegion;
    delete callableSpan;
    delete maximumSpan;
    delete extendedSpan;
    delete idealSpan;
    delete nonVariantFlanks->first;
    delete nonVariantFlanks->second;
    delete nonVariantFlanks;
}

bool AssemblyRegionTrimmer_Result::isVariationPresent() {
    return !callableEvents.empty();
}

AssemblyRegion *AssemblyRegionTrimmer_Result::getCallableRegion() {
    if(callableRegion == nullptr && extendedSpan == nullptr) {
        callableRegion = emitReferenceConfidence ? originalRegion
    }
}
