//
// Created by 梦想家xixi on 2021/10/14.
//

#include "AssemblyRegion.h"
#include "IntervalUtils.h"

AssemblyRegion::AssemblyRegion(SimpleInterval& activeRegionLoc,
                               std::vector<ActivityProfileState> supportingStates, const bool isActivate,
                               const int extension) : activeRegionLoc(activeRegionLoc), supportingStates(std::move(supportingStates)), isActivate(isActivate), extension(extension){
    std::string contig = activeRegionLoc.getContig();
    SimpleInterval* simpleInterval = trimIntervalToContig(contig, activeRegionLoc.getStart() - extension, activeRegionLoc.getEnd() + extension);
    extendedLoc = *simpleInterval;
    spanIncludingReads = extendedLoc;
    delete simpleInterval;
}

SimpleInterval *AssemblyRegion::trimIntervalToContig(std::string& contig, const int start, const int stop) {
//    const int contigLength = header.getSequence(contig).getSequenceLength();
    return IntervalUtils::trimIntervalToContig(contig, start, stop, 10);
}

void AssemblyRegion::checkStates(SimpleInterval &activeRegion) {
    if(!supportingStates.empty()) {
        Mutect2Utils::validateArg(supportingStates.size() == activeRegionLoc.size(), "Supporting states wasn't empty but it doesn't have exactly one state per bp in the active region.");
        std::vector<ActivityProfileState>::iterator pr;
        for(pr = supportingStates.begin(); pr != supportingStates.end() - 1; pr++) {
            Mutect2Utils::validateArg((pr+1)->getLoc().getStart() == pr->getLoc().getStart() + 1 &&
                                              (pr+1)->getLoc().getContig() == pr->getLoc().getContig(),
                                          "Supporting state has an invalid sequence");
        }
    }
}

AssemblyRegion::~AssemblyRegion() = default;

