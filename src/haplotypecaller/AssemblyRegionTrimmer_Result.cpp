//
// Created by 梦想家xixi on 2021/12/14.
//

#include "AssemblyRegionTrimmer_Result.h"

#include <utility>

AssemblyRegionTrimmer_Result::AssemblyRegionTrimmer_Result(bool emitReferenceConfidence, bool needsTrimming,
                                                           const std::shared_ptr<AssemblyRegion>& originalRegion, int padding, int extension,
                                                           std::vector<std::shared_ptr<VariantContext>> * overlappingEvents,
                                                           std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>  nonVariantFlanks,
                                                           const std::shared_ptr<SimpleInterval>& extendedSpan, const std::shared_ptr<SimpleInterval>& idealSpan,
                                                           const std::shared_ptr<SimpleInterval>& maximumSpan, const std::shared_ptr<SimpleInterval>& callableSpan) : emitReferenceConfidence(emitReferenceConfidence), needsTrimming(needsTrimming), callableEvents(*overlappingEvents),padding(padding),
                                                                                                                        usableExtension(extension), nonVariantFlanks(std::move(nonVariantFlanks)), extendedSpan(extendedSpan), idealSpan(idealSpan), maximumSpan(maximumSpan), callableSpan(callableSpan),
                                                                                                                        originalRegion(originalRegion)
                                                          {
    Mutect2Utils::validateArg(extendedSpan == nullptr || callableSpan == nullptr || extendedSpan->contains(callableSpan), "the extended callable span must include the callable span");

}

std::shared_ptr<AssemblyRegionTrimmer_Result>
AssemblyRegionTrimmer_Result::noVariation(bool emitReferenceConfidence, std::shared_ptr<AssemblyRegion> targetRegion, int padding,
                                          int usableExtension) {
    std::vector<std::shared_ptr<VariantContext>> events;
    std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>> nonVariantFlanks = std::make_shared<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>(targetRegion->getSpan(), nullptr);
    std::shared_ptr<AssemblyRegionTrimmer_Result> result(new AssemblyRegionTrimmer_Result(emitReferenceConfidence, false, targetRegion, padding, usableExtension,
                                                                            &events, nonVariantFlanks,
                                                                            nullptr, nullptr, nullptr, nullptr));
    result->leftFlankRegion = targetRegion;
    return result;
}

std::shared_ptr<AssemblyRegionTrimmer_Result>
AssemblyRegionTrimmer_Result::noTrimming(bool emitReferenceConfidence, const std::shared_ptr<AssemblyRegion>& targetRegion, int padding,
                                         int usableExtension, std::vector<std::shared_ptr<VariantContext>> *events) {
    std::shared_ptr<SimpleInterval> targetRegionLoc = targetRegion->getSpan();
    std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>> nonVariantFlanks = std::make_shared<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>(nullptr, nullptr);
    std::shared_ptr<AssemblyRegionTrimmer_Result> result(new AssemblyRegionTrimmer_Result(emitReferenceConfidence, false, targetRegion, padding, usableExtension, events, nonVariantFlanks, targetRegionLoc, targetRegionLoc, targetRegionLoc, targetRegionLoc));
    result->callableRegion = targetRegion;
    return result;
}

AssemblyRegionTrimmer_Result::~AssemblyRegionTrimmer_Result() {
}

bool AssemblyRegionTrimmer_Result::isVariationPresent() {
    return !callableEvents.empty();
}

std::shared_ptr<AssemblyRegion> AssemblyRegionTrimmer_Result::getCallableRegion() {
    if(callableRegion == nullptr && extendedSpan == nullptr) {
        callableRegion = emitReferenceConfidence ? originalRegion->trim(callableSpan, extendedSpan) : originalRegion->trim(extendedSpan, extendedSpan);
    }
    return callableRegion;
}

bool AssemblyRegionTrimmer_Result::getNeedsTrimming() {
    return needsTrimming;
}
