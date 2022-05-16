//
// Created by 梦想家xixi on 2021/10/14.
//

#include "AssemblyRegion.h"
#include <cassert>
#include <utility>
#include "IntervalUtils.h"
#include "clipping/ReadClipper.h"

AssemblyRegion::AssemblyRegion(SimpleInterval const &activeRegionLoc,
                               std::vector<std::shared_ptr<ActivityProfileState>> supportingStates, const bool isActive,
                               const int extension, SAMFileHeader *header)
		: activeRegionLoc(std::make_shared<SimpleInterval>(activeRegionLoc)),
		  supportingStates(std::move(supportingStates)), isActive(isActive), extension(extension), header(header) {

	std::string contig = activeRegionLoc.getContig();
	extendedLoc = trimIntervalToContig(contig, activeRegionLoc.getStart() - extension,
	                                   activeRegionLoc.getEnd() + extension);
	spanIncludingReads = extendedLoc;
	//checkStates(this->activeRegionLoc);
}

AssemblyRegion::AssemblyRegion(SimpleInterval const &activeRegionLoc, const int extension)
		: activeRegionLoc(std::make_shared<SimpleInterval>(activeRegionLoc)), isActive(true), extension(extension) {
	std::string contig = activeRegionLoc.getContig();
	extendedLoc = trimIntervalToContig(contig, activeRegionLoc.getStart() - extension,
	                                   activeRegionLoc.getEnd() + extension);
	spanIncludingReads = extendedLoc;
	//checkStates(this->activeRegionLoc);
}

std::shared_ptr<SimpleInterval>
AssemblyRegion::trimIntervalToContig(std::string &contig, const int start, const int stop) {
	const int contigLength = header->getSequenceDictionary().getSequence(contig).getSequenceLength();
	return IntervalUtils::trimIntervalToContig(contig, start, stop, contigLength);
}

void AssemblyRegion::checkStates(SimpleInterval &activeRegion) {
	if (!supportingStates.empty()) {
		Mutect2Utils::validateArg(supportingStates.size() == activeRegionLoc->size(),
		                          "Supporting states wasn't empty but it doesn't have exactly one state per bp in the active region.");
		std::vector<std::shared_ptr<ActivityProfileState>>::iterator pr;
		for (pr = supportingStates.begin(); pr != supportingStates.end() - 1; pr++) {
			Mutect2Utils::validateArg((*(pr + 1))->getLoc().getStart() == (*pr)->getLoc().getStart() + 1 &&
			                          (*(pr + 1))->getLoc().getContig() == (*pr)->getLoc().getContig(),
			                          "Supporting state has an invalid sequence");
		}
	}
}

AssemblyRegion::~AssemblyRegion()
{
/*    for(int i=0; i<reads.size(); i++)
    {
        std::cout << i << " " << reads[i].use_count() << std::endl;
    }*/
}

std::ostream &operator<<(std::ostream &os, AssemblyRegion &assemblyRegion) {
	os << "AssemblyRegion " << *assemblyRegion.activeRegionLoc << "active:   " << assemblyRegion.isActive << std::endl;
	return os;
}

void AssemblyRegion::setIsActive(const bool value) { isActive = value; }

bool AssemblyRegion::equalsIgnoreReads(const AssemblyRegion &other) {
	if (!(activeRegionLoc == other.activeRegionLoc) || isActive != other.getIsActive() ||
	    extension != other.getExtension())
		return false;
	return extendedLoc == activeRegionLoc;
}

void AssemblyRegion::setFinalized(bool value) { hasBeenFinalized = value; }

SAMFileHeader *AssemblyRegion::getHeader() {
	return header;
}

std::vector<std::shared_ptr<SAMRecord>> &AssemblyRegion::getReads() {
	return reads;
}

void AssemblyRegion::setRead(std::vector<std::shared_ptr<SAMRecord>> &reads) {
	this->reads = reads;
}

std::shared_ptr<AssemblyRegion>
AssemblyRegion::trim(const std::shared_ptr<SimpleInterval> &span, const std::shared_ptr<SimpleInterval> &extendedSpan) {
	Mutect2Utils::validateArg(span.get(), "Active region extent cannot be null");
	Mutect2Utils::validateArg(extendedSpan.get(), "Active region extended span cannot be null");
	Mutect2Utils::validateArg(extendedSpan->contains(span),
	                          "The requested extended span must fully contain the requested span");

	std::shared_ptr<SimpleInterval> subActive = getSpan()->intersect(span);
	int requiredOnRight = std::max(extendedSpan->getEnd() - subActive->getEnd(), 0);
	int requiredOnLeft = std::max(subActive->getStart() - extendedSpan->getStart(), 0);
	int requiredExtension = std::min(std::max(requiredOnLeft, requiredOnRight), getExtension());

	std::shared_ptr<AssemblyRegion> result = std::make_shared<AssemblyRegion>(*subActive,
	                                                                          std::vector<std::shared_ptr<ActivityProfileState>>(),
	                                                                          isActive, requiredExtension, header);
	std::vector<std::shared_ptr<SAMRecord>> myReads = getReads();
	std::shared_ptr<SimpleInterval> resultExtendedLoc = result->getExtendedSpan();
	int resultExtendedLocStart = resultExtendedLoc->getStart();
	int resultExtendedLocStop = resultExtendedLoc->getEnd();

	std::vector<std::shared_ptr<SAMRecord>> trimmedReads;
	for (const std::shared_ptr<SAMRecord> &read: myReads) {
		std::shared_ptr<SAMRecord> clippedRead = ReadClipper::hardClipToRegion(read, resultExtendedLocStart,
		                                                                       resultExtendedLocStop);
		if (result->readOverlapsRegion(clippedRead) && !clippedRead->isEmpty()) {
			trimmedReads.emplace_back(clippedRead);
		}
	}
	//TODO:sort
	result->clearReads();
	result->addAll(trimmedReads);
	return result;
}

bool AssemblyRegion::readOverlapsRegion(std::shared_ptr<SAMRecord> &read) {
	if (read->isEmpty() || read->getStart() > read->getEnd()) {
		return false;
	}
	std::shared_ptr<SimpleInterval> readLoc = std::make_shared<SimpleInterval>(read->getContig(), read->getStart(),
	                                                                           read->getEnd());
	return readLoc->overlaps(extendedLoc);
}

void AssemblyRegion::clearReads() {
	spanIncludingReads = extendedLoc;
	reads.clear();
}

void AssemblyRegion::addAll(std::vector<std::shared_ptr<SAMRecord>> &readsToAdd) {
	for (std::shared_ptr<SAMRecord> &samRecord: readsToAdd) {
		reads.emplace_back(samRecord);
	}
}

std::shared_ptr<uint8_t[]> AssemblyRegion::getAssemblyRegionReference(ReferenceCache *cache, int padding, int &length) {
	return getReference(cache, padding, extendedLoc, length);
}

std::shared_ptr<uint8_t[]> AssemblyRegion::getReference(ReferenceCache *referenceReader, int padding,
                                                        const std::shared_ptr<SimpleInterval> &genomeLoc, int &length) {
	Mutect2Utils::validateArg(referenceReader, "referenceReader cannot be null");
	Mutect2Utils::validateArg(padding >= 0, "padding must be a positive integer but got");
	Mutect2Utils::validateArg(genomeLoc->size() > 0, "GenomeLoc must have size > 0 but got ");
	std::string contig = genomeLoc->getContig();
	return referenceReader->getSubsequenceAt(header->getSequenceDictionary().getSequenceIndex(contig),
	                                         std::max(0, genomeLoc->getStart() - padding),
	                                         std::min(header->getSequenceDictionary().getSequence(
			                                                  contig).getSequenceLength(),
	                                                  genomeLoc->getEnd() + padding), length);
}

void AssemblyRegion::removeAll(const std::vector<std::shared_ptr<SAMRecord>> &readsToRemove) {
	for (const auto &read: readsToRemove) {
		reads.erase(std::find(reads.begin(), reads.end(), read));
	}
}

void AssemblyRegion::add(std::shared_ptr<SAMRecord> &read) {
    assert(read != nullptr);
    assert(readOverlapsRegion(read));

    spanIncludingReads = spanIncludingReads->mergeWithContiguous(read->getLoc());

    if(!reads.empty())
    {
        std::shared_ptr<SAMRecord> & lastRead = reads.back();
        assert(lastRead->getContig() == read->getContig());
        //assert(read->getStart() >= lastRead->getStart());  // TODO: Is this line necessary?
    }
    reads.emplace_back(std::make_shared<SAMRecord>(*read));
}
