//
// Created by 梦想家xixi on 2021/12/14.
//

#include "AssemblyRegionTrimmer_Result.h"

AssemblyRegionTrimmer_Result::AssemblyRegionTrimmer_Result(bool emitReferenceConfidence, bool needsTrimming,
                                                           AssemblyRegion *originalRegion, int padding, int extension,
                                                           std::vector<VariantContext *> * overlappingEvents,
                                                           std::pair<SimpleInterval *, SimpleInterval *> * nonVariantFlanks,
                                                           SimpleInterval *extendedSpan, SimpleInterval *idealSpan,
                                                           SimpleInterval *maximumSpan, SimpleInterval *callableSpan) : emitReferenceConfidence(emitReferenceConfidence), needsTrimming(needsTrimming), originalRegion(originalRegion),
                                                           nonVariantFlanks(nonVariantFlanks), padding(padding), usableExtension(extension), callableEvents(overlappingEvents), idealSpan(idealSpan), maximumSpan(maximumSpan), extendedSpan(extendedSpan){
    Mutect2Utils::validateArg(extendedSpan == nullptr || callableSpan == nullptr || extendedSpan->contains(callableSpan), "the extended callable span must include the callable span");
}
