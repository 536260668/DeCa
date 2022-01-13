//
// Created by 梦想家xixi on 2021/10/14.
//

#include "assert.h"
#include "AssemblyRegion.h"
#include "IntervalUtils.h"
#include "clipping/ReadClipper.h"

AssemblyRegion::AssemblyRegion(SimpleInterval const &activeRegionLoc,
                               std::vector<ActivityProfileState> supportingStates, const bool isActive,
                               const int extension, sam_hdr_t * header) : activeRegionLoc(activeRegionLoc), supportingStates(std::move(supportingStates)), isActive(isActive),
                               extension(extension), hdr(header){

    std::string contig = activeRegionLoc.getContig();
    SimpleInterval* simpleInterval = trimIntervalToContig(contig, activeRegionLoc.getStart() - extension, activeRegionLoc.getEnd() + extension);
    assert(simpleInterval != nullptr);
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
    //const int contigLength = sam_hdr_tid2len(hdr, sam_hdr_name2tid(hdr, contig.c_str()));
    int contigLength = 1000000;
    return IntervalUtils::trimIntervalToContig(contig, start, stop, contigLength);
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
    if(!(activeRegionLoc == other.activeRegionLoc) || isActive != other.getIsActive() || extension != other.getExtension())
        return false;
    return extendedLoc == activeRegionLoc;
}
void AssemblyRegion::setFinalized(bool value) {hasBeenFinalized = value;}

sam_hdr_t * AssemblyRegion::getHeader()
{
    return hdr;
}

std::vector<std::shared_ptr<SAMRecord>> & AssemblyRegion::getReads(){
    return reads;
}

void AssemblyRegion::setRead(std::vector<std::shared_ptr<SAMRecord>> &reads) {
    this->reads = reads;
}

AssemblyRegion *AssemblyRegion::trim(SimpleInterval *span, SimpleInterval *extendedSpan) {
    Mutect2Utils::validateArg(span, "Active region extent cannot be null");
    Mutect2Utils::validateArg(extendedSpan, "Active region extended span cannot be null");
    Mutect2Utils::validateArg(extendedSpan->contains(span), "The requested extended span must fully contain the requested span");

    SimpleInterval* subActive = getSpan().intersect(span);
    int requiredOnRight = std::max(extendedSpan->getEnd() - subActive->getEnd(), 0);
    int requiredOnLeft = std::max(subActive->getStart() - extendedSpan->getStart(), 0);
    int requiredExtension = std::min(std::max(requiredOnLeft, requiredOnRight), getExtension());

    AssemblyRegion* result = new AssemblyRegion(subActive, std::vector<ActivityProfileState>(), isActive, requiredExtension, hdr);
    std::vector<std::shared_ptr<SAMRecord>>  myReads = getReads();
    SimpleInterval resultExtendedLoc = result->getExtendedSpan();
    int resultExtendedLocStart = resultExtendedLoc.getStart();
    int resultExtendedLocStop = resultExtendedLoc.getEnd();

    std::vector<std::shared_ptr<SAMRecord>> trimmedReads;
    for(std::shared_ptr<SAMRecord> read : myReads) {
        std::shared_ptr<SAMRecord> clippedRead = ReadClipper::hardClipToRegion(read, resultExtendedLocStart, resultExtendedLocStop);
        if(result->readOverlapsRegion(clippedRead) && !clippedRead->isEmpty()) {
            trimmedReads.emplace_back(clippedRead);
        }
    }
    //TODO:sort
    result->clearReads();
    result->addAll(trimmedReads);
    return result;
}

bool AssemblyRegion::readOverlapsRegion(std::shared_ptr<SAMRecord> & read) {
    if(read->isEmpty() || read->getStart() > read->getEnd()) {
        return false;
    }
    SimpleInterval readLoc(read->getContig(), read->getStart(), read->getEnd());
    return readLoc.overlaps(&extendedLoc);
}

void AssemblyRegion::clearReads() {
    spanIncludingReads = extendedLoc;
    reads.clear();
}

void AssemblyRegion::addAll(std::vector<std::shared_ptr<SAMRecord>> &readsToAdd) {
    for(std::shared_ptr<SAMRecord> & samRecord : readsToAdd) {
        reads.emplace_back(samRecord);
    }
}
