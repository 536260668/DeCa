//
// Created by 梦想家xixi on 2021/12/8.
//

#include <limits>
#include "AssemblyBasedCallerUtils.h"
#include "haplotypecaller/ReferenceConfidenceModel.h"
#include "clipping/ReadClipper.h"
#include "read/ReadUtils.h"
#include "QualityUtils.h"
#include "utils/fragments/FragmentCollection.h"
#include "utils/fragments/FragmentUtils.h"
#include "AlignmentUtils.h"

std::shared_ptr<Haplotype>
AssemblyBasedCallerUtils::createReferenceHaplotype(const std::shared_ptr<AssemblyRegion> &region,
                                                   const std::shared_ptr<SimpleInterval> &referencePadding,
                                                   ReferenceCache &cache) {
	int length = 0;
	std::shared_ptr<uint8_t[]> tmp = region->getAssemblyRegionReference(&cache, 0, length);
	return ReferenceConfidenceModel::createReferenceHaplotype(region, tmp, length, referencePadding);
}

std::shared_ptr<AssemblyResultSet>
AssemblyBasedCallerUtils::assembleReads(const std::shared_ptr<AssemblyRegion> &region,
                                        M2ArgumentCollection &argumentCollection,
                                        SAMFileHeader *header, ReferenceCache &cache,
                                        ReadThreadingAssembler &assemblyEngine) {
	finalizeRegion(region, false, false, 9, header, false);
	int refLength = 0;
	std::shared_ptr<uint8_t[]> fullReferenceWithPadding = region->getAssemblyRegionReference(&cache,
	                                                                                         REFERENCE_PADDING_FOR_ASSEMBLY,
	                                                                                         refLength);
	const std::shared_ptr<SimpleInterval> paddedReferenceLoc = getPaddedReferenceLoc(region,
	                                                                                 REFERENCE_PADDING_FOR_ASSEMBLY,
	                                                                                 header);
	std::shared_ptr<Haplotype> refHaplotype = createReferenceHaplotype(region, paddedReferenceLoc, cache);
	std::shared_ptr<AssemblyResultSet> assemblyResultSet = assemblyEngine.runLocalAssembly(region, refHaplotype,
	                                                                                       fullReferenceWithPadding,
	                                                                                       refLength,
	                                                                                       paddedReferenceLoc,
	                                                                                       nullptr);
	return assemblyResultSet;
}

std::shared_ptr<SimpleInterval>
AssemblyBasedCallerUtils::getPaddedReferenceLoc(const std::shared_ptr<AssemblyRegion> &region, int referencePadding,
                                                SAMFileHeader *header) {
	int padLeft = std::max(region->getExtendedSpan()->getStart() - referencePadding, 0);
	int padRight = std::min(region->getExtendedSpan()->getEnd() + referencePadding,
	                        header->getSequenceDictionary().getSequence(
			                        region->getExtendedSpan()->getContig()).getSequenceLength() - 1);
	return std::make_shared<SimpleInterval>(region->getExtendedSpan()->getContig(), padLeft, padRight);
}

void
AssemblyBasedCallerUtils::finalizeRegion(const std::shared_ptr<AssemblyRegion> &region, bool errorCorrectReads,
                                         bool dontUseSoftClippedBases,
                                         uint8_t minTailQuality, SAMFileHeader *header,
                                         bool correctOverlappingBaseQualities) {
	if (region->isFinalized())
		return;
	std::vector<std::shared_ptr<SAMRecord>> readsToUse;
	for (const std::shared_ptr<SAMRecord> &myRead: region->getReads()) {
		uint8_t minTailQualityToUse = errorCorrectReads ? 6 : minTailQuality;
		std::shared_ptr<SAMRecord> clippedRead = ReadClipper::hardClipLowQualEnds(myRead, minTailQuality);
		clippedRead = dontUseSoftClippedBases || !ReadUtils::hasWellDefinedFragmentSize(clippedRead) ?
		              ReadClipper::hardClipSoftClippedBases(clippedRead) : ReadClipper::revertSoftClippedBases(
						clippedRead);

		clippedRead = clippedRead->isUnmapped() ? clippedRead : ReadClipper::hardClipAdaptorSequence(clippedRead);
		if (!clippedRead->isEmpty() && clippedRead->getCigar()->getReadLength() > 0) {
			clippedRead = ReadClipper::hardClipToRegion(clippedRead, region->getExtendedSpan()->getStart(),
			                                            region->getExtendedSpan()->getEnd());
			if (region->readOverlapsRegion(clippedRead) && clippedRead->getLength() > 0) {
				readsToUse.emplace_back(
						(clippedRead == myRead) ? std::shared_ptr<SAMRecord>(new SAMRecord(*clippedRead))
						                        : clippedRead);
			}
		}
	}
	region->clearReads();
	region->addAll(readsToUse);
	region->setFinalized(true);
}

std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>
AssemblyBasedCallerUtils::splitReadsBySample(const std::vector<std::shared_ptr<SAMRecord>> &reads) {
	std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> res = std::make_shared<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>();
	std::string normal = "normal";
	std::string tumor = "case";     // TODO: make it a parameter, not a constant string
	res->insert({normal, std::vector<std::shared_ptr<SAMRecord>>()});
	res->insert({tumor, std::vector<std::shared_ptr<SAMRecord>>()});
	std::vector<std::shared_ptr<SAMRecord>> &normalReads = res->at(normal);
	std::vector<std::shared_ptr<SAMRecord>> &tumorReads = res->at(tumor);
	for (const std::shared_ptr<SAMRecord> &read: reads) {
		if (read->getGroup() == 0) {
			normalReads.emplace_back(read);
		} else {
			tumorReads.emplace_back(read);
		}
	}
	return res;
}

PairHMMLikelihoodCalculationEngine *
AssemblyBasedCallerUtils::createLikelihoodCalculationEngine(LikelihoodEngineArgumentCollection& likelihoodArgs)
{
    double log10GlobalReadMismappingRate = likelihoodArgs.phredScaledGlobalReadMismappingRate < 0 ? (-1) * std::numeric_limits<double>::infinity()
            : QualityUtils::qualToErrorProbLog10(likelihoodArgs.phredScaledGlobalReadMismappingRate);
    return new PairHMMLikelihoodCalculationEngine((char)likelihoodArgs.gcpHMM, likelihoodArgs.pairHMMNativeArgs, log10GlobalReadMismappingRate, likelihoodArgs.pcrErrorModel, likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD);
}

void AssemblyBasedCallerUtils::cleanOverlappingReadPairs(vector<shared_ptr<SAMRecord>> &reads, string sample,
                                                         bool setConflictingToZero, int halfOfPcrSnvQual,
                                                         int halfOfPcrIndelQual) {
    auto MappedReads = splitReadsBySample(reads);
    for(auto iter = MappedReads->begin(); iter != MappedReads->end(); iter++)
    {
        FragmentCollection<SAMRecord> * fragmentCollection = FragmentCollection<SAMRecord>::create(iter->second);
        for(std::pair<std::shared_ptr<SAMRecord>, std::shared_ptr<SAMRecord>>& overlappingPair : fragmentCollection->getOverlappingPairs())
        {
            FragmentUtils::adjustQualsOfOverlappingPairedFragments(overlappingPair, setConflictingToZero, halfOfPcrSnvQual, halfOfPcrIndelQual);
        }
        delete fragmentCollection;
    }
}

std::shared_ptr<AssemblyRegion> AssemblyBasedCallerUtils::assemblyRegionWithWellMappedReads(const std::shared_ptr<AssemblyRegion>& originalAssemblyRegion, int minMappingQuality, SAMFileHeader * header)
{
    auto result = make_shared<AssemblyRegion>(*originalAssemblyRegion->getSpan(), originalAssemblyRegion->getSupportingStates(), originalAssemblyRegion->getIsActive(), originalAssemblyRegion->getExtension(), header);
    for(auto & read : originalAssemblyRegion->getReads())
    {
        if(read->getMappingQuality() >= minMappingQuality)
            result->add(read);
    }
    return result;
}

shared_ptr<unordered_map<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>> AssemblyBasedCallerUtils::realignReadsToTheirBestHaplotype(
        AlleleLikelihoods<SAMRecord, Haplotype> &originalReadLikelihoods, shared_ptr<Haplotype> &refHaplotype, shared_ptr<SimpleInterval>& paddedReferenceLoc,
        SmithWatermanAligner *aligner) {
    auto bestAlleles = originalReadLikelihoods.bestAllelesBreakingTies(&AssemblyBasedCallerUtils::HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY);
    auto result = std::make_shared<unordered_map<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>>(bestAlleles->size());

    for(auto & bestAllele : *bestAlleles)
    {
        auto & originalRead = bestAllele->evidence;
        auto & bestHaplotype = bestAllele->allele;
        bool isInformative = bestAllele->isInformative();
        shared_ptr<SAMRecord> realignedRead = AlignmentUtils::createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc->getStart(), isInformative, aligner);
        result->insert(pair<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>(originalRead, realignedRead));
    }
    return result;
}

double AssemblyBasedCallerUtils::HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY(shared_ptr<Haplotype> h) {
    assert(h != nullptr);
    auto cigar = h->getCigar();
    int referenceTerm = h->getIsReference() ? 1 : 0;
    int cigarTerm = cigar == nullptr ? 0 : 1-cigar->numCigarElements();

    return (double)(referenceTerm + cigarTerm);
}