//
// Created by 梦想家xixi on 2021/12/14.
//

#include "AssemblyRegionTrimmer.h"

AssemblyRegionTrimmer::AssemblyRegionTrimmer(ReadThreadingAssemblerArgumentCollection *assemblyArgs,
                                             SAMSequenceDictionary *sequenceDictionary, bool isGGA,
                                             bool emitReferenceConfidence) : assemblyArgs(assemblyArgs), sequenceDictionary(sequenceDictionary),
                                             emitReferenceConfidence(emitReferenceConfidence){
    Mutect2Utils::validateArg(assemblyArgs != nullptr, "null is not allowed there");
    checkUserArguments();
    usableExtension = isGGA ? assemblyArgs->ggaExtension : assemblyArgs->discoverExtension;
}

void AssemblyRegionTrimmer::checkUserArguments() {
    if(assemblyArgs->snpPadding < 0) {
        throw std::invalid_argument("paddingAroundSNPs");
    }
    if(assemblyArgs->indelPadding < 0) {
        throw std::invalid_argument("paddingAroundIndels");
    }
    if(assemblyArgs->discoverExtension < 0) {
        throw std::invalid_argument("maxDiscARExtension");
    }
    if(assemblyArgs->ggaExtension < 0) {
        throw std::invalid_argument("maxGGAAREExtension");
    }
}

std::shared_ptr<AssemblyRegionTrimmer_Result> AssemblyRegionTrimmer::trim(const std::shared_ptr<AssemblyRegion>& originalRegion,
                                                          std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & allVariantsWithinExtendedRegion) {
    if(allVariantsWithinExtendedRegion.empty()) {
        return AssemblyRegionTrimmer_Result::noVariation(emitReferenceConfidence, originalRegion, assemblyArgs->snpPadding, usableExtension);
    }
    std::vector<std::shared_ptr<VariantContext>> withinActiveRegion;
    SimpleInterval & originalRegionRange = originalRegion->getSpan();
    bool foundNonSnp = false;
    SimpleInterval variantSpan;
    bool flag = false;
    for(std::shared_ptr<VariantContext> vc : allVariantsWithinExtendedRegion) {
        SimpleInterval vcLoc(vc->getContig(), vc->getStart(), vc->getEnd());
        if(originalRegionRange.overlaps(&vcLoc)) {
            foundNonSnp = foundNonSnp || !vc->isSNP();
            variantSpan = !flag ? vcLoc : variantSpan.spanWith(&vcLoc);
            if(!flag){
                flag = true;
            }
            withinActiveRegion.emplace_back(vc);
        }
    }
    int padding = foundNonSnp ? assemblyArgs->indelPadding : assemblyArgs->snpPadding;
    if(variantSpan == SimpleInterval()) {
        return AssemblyRegionTrimmer_Result::noVariation(emitReferenceConfidence, originalRegion, padding, usableExtension);
    }

    if(assemblyArgs->dontTrimActiveRegions) {
        return AssemblyRegionTrimmer_Result::noTrimming(emitReferenceConfidence, originalRegion, padding, usableExtension, &withinActiveRegion);
    }
    SimpleInterval* maximumSpan = originalRegionRange.expandWithinContig(usableExtension, sequenceDictionary);
    SimpleInterval* idealSpan = variantSpan.expandWithinContig(padding, sequenceDictionary);
    SimpleInterval* finalSpan = maximumSpan->intersect(idealSpan)->mergeWithContiguous(&variantSpan);
    SimpleInterval* callableSpan = emitReferenceConfidence ? variantSpan.intersect(&originalRegionRange) : &variantSpan;
    std::pair<SimpleInterval *, SimpleInterval *> * nonVariantRegions = nonVariantTargetRegions(originalRegion, callableSpan);
    std::shared_ptr<AssemblyRegionTrimmer_Result> ret(new AssemblyRegionTrimmer_Result(emitReferenceConfidence, true, originalRegion, padding, usableExtension, &withinActiveRegion, nonVariantRegions, finalSpan, idealSpan, maximumSpan, &variantSpan));
    delete maximumSpan;
    delete idealSpan;
    delete finalSpan;
    if(emitReferenceConfidence)
        delete callableSpan;
    delete nonVariantRegions->first;
    delete nonVariantRegions->second;
    delete nonVariantRegions;
    return ret;
}

std::pair<SimpleInterval *, SimpleInterval *> *
AssemblyRegionTrimmer::nonVariantTargetRegions(std::shared_ptr<AssemblyRegion> targetRegion, SimpleInterval *variantSpan) {
    SimpleInterval targetRegionRange = targetRegion->getSpan();
    int finalStart = variantSpan->getStart();
    int finalStop = variantSpan->getEnd();
    int targetStart = targetRegionRange.getStart();
    int targetStop = targetRegionRange.getEnd();
    bool preTrimmingRequired = targetStart < finalStart;
    bool postTrimmingRequired = targetStop > finalStop;
    if(preTrimmingRequired) {
        std::string contig = targetRegionRange.getContig();
        return postTrimmingRequired ? new std::pair<SimpleInterval*, SimpleInterval*>(new SimpleInterval(contig, targetStart, finalStart-1), new SimpleInterval(contig, finalStop+1, targetStop)) :
                new std::pair<SimpleInterval*, SimpleInterval*>(new SimpleInterval(contig, targetStart, finalStart-1), nullptr);
    } else if (postTrimmingRequired) {
        return new std::pair<SimpleInterval*, SimpleInterval*>(nullptr, new SimpleInterval(targetRegionRange.getContig(), finalStop+1, targetStop));
    } else {
        return new std::pair<SimpleInterval*, SimpleInterval*>(nullptr, nullptr);
    }
}
