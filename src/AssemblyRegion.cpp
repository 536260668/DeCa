//
// Created by 梦想家xixi on 2021/10/14.
//

#include "assert.h"
#include "AssemblyRegion.h"
#include "IntervalUtils.h"
#include "clipping/ReadClipper.h"

AssemblyRegion::AssemblyRegion(SimpleInterval const &activeRegionLoc,
                               std::vector<ActivityProfileState> supportingStates, const bool isActive,
                               const int extension, SAMFileHeader * header) : activeRegionLoc(activeRegionLoc), supportingStates(std::move(supportingStates)), isActive(isActive),
                               extension(extension), header(header){

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
    const int contigLength = header->getSequenceDictionary().getSequence(contig).getSequenceLength();
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

SAMFileHeader * AssemblyRegion::getHeader()
{
    return header;
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

    AssemblyRegion* result = new AssemblyRegion(subActive, std::vector<ActivityProfileState>(), isActive, requiredExtension, header);
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

uint8_t *AssemblyRegion::getAssemblyRegionReference(ReferenceCache *cache, int padding, int & length) {
    return getReference(cache, padding, extendedLoc, length);
}

uint8_t *AssemblyRegion::getReference(ReferenceCache *referenceReader, int padding, SimpleInterval &genomeLoc, int & length) {
    Mutect2Utils::validateArg(referenceReader, "referenceReader cannot be null");
    Mutect2Utils::validateArg(padding >= 0, "padding must be a positive integer but got");
    Mutect2Utils::validateArg(genomeLoc.size() > 0, "GenomeLoc must have size > 0 but got ");
    std::string contig = genomeLoc.getContig();
    return referenceReader->getSubsequenceAt(header->getSequenceDictionary().getSequenceIndex(contig), std::max(0, genomeLoc.getStart() - padding),
                                             std::min(header->getSequenceDictionary().getSequence(contig).getSequenceLength(), genomeLoc.getEnd() + padding), length);
}
