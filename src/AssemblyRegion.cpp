//
// Created by 梦想家xixi on 2021/10/14.
//

#include "AssemblyRegion.h"
#include "IntervalUtils.h"

AssemblyRegion::AssemblyRegion(SimpleInterval const &activeRegionLoc,
                               std::vector<ActivityProfileState> supportingStates, const bool isActive,
                               const int extension) : activeRegionLoc(activeRegionLoc), supportingStates(std::move(supportingStates)), isActive(isActive), extension(extension){
    std::string contig = activeRegionLoc.getContig();
    SimpleInterval* simpleInterval = trimIntervalToContig(contig, activeRegionLoc.getStart() - extension, activeRegionLoc.getEnd() + extension);
    extendedLoc = *simpleInterval;
    spanIncludingReads = extendedLoc;
    delete simpleInterval;
    checkStates(this->activeRegionLoc);
}

AssemblyRegion::AssemblyRegion(SimpleInterval const &activeRegionLoc, const int extension) : activeRegionLoc(activeRegionLoc),  isActive(
        true), extension(extension){
    std::string contig = activeRegionLoc.getContig();
    SimpleInterval* simpleInterval = trimIntervalToContig(contig, activeRegionLoc.getStart() - extension, activeRegionLoc.getEnd() + extension);
    extendedLoc = *simpleInterval;
    spanIncludingReads = extendedLoc;
    delete simpleInterval;
    checkStates(this->activeRegionLoc);
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

std::ostream & operator<<(std::ostream & os, AssemblyRegion & assemblyRegion) {
    os << "AssemblyRegion " << assemblyRegion.activeRegionLoc << "active:   " << assemblyRegion.isActive << std::endl;
    return os;
}

void AssemblyRegion::setIsActive(const bool value) {isActive = value;}

bool AssemblyRegion::equalsIgnoreReads(const AssemblyRegion &other) {
    if(!(activeRegionLoc == other.getSpan()) || isActive != other.getIsActive() || extension != other.getExtension())
        return false;
    return extendedLoc == other.getExtendedSpan();
}
void AssemblyRegion::setFinalized(bool value) {hasBeenFinalized = value;}

