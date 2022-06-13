//
// Created by lhh on 10/19/21.
//

#include <cmath>
#include <utility>
#include <cassert>
#include "Mutect2Engine.h"
#include "cigar/Cigar.h"
#include "utils/PeUtils.h"
#include "samtools/SAMRecord.h"
#include "MathUtils.h"
#include "QualityUtils.h"
#include "NaturalLogUtils.h"
#include "AssemblyResultSet.h"
#include "haplotypecaller/AssemblyBasedCallerUtils.h"

Mutect2Engine::Mutect2Engine(M2ArgumentCollection & MTAC, SAMFileHeader* samFileHeader, const std::string &modelPath, VariantAnnotatorEngine& annotatorEngine):MTAC(MTAC), minCallableDepth(MTAC.callableDepth),
                                                            normalSample(MTAC.normalSample) ,callableSites(0), refCache(nullptr) ,header(samFileHeader),
                                                                                                    assemblyEngine(0, 1, 128, false, false, {10, 25}),
                                                                                                    likelihoodCalculationEngine(AssemblyBasedCallerUtils::createLikelihoodCalculationEngine(MTAC.likelihoodArgs)),
                                                                                                    trimmer(&assemblerArgs, &header->getSequenceDictionary(), false,
                                                                                                            false), aligner(SmithWatermanAligner::getAligner(SmithWatermanAligner::FASTEST_AVAILABLE)),
                                                                                                            genotypingEngine(MTAC, MTAC.normalSample, annotatorEngine)
{
    std::vector<SAMReadGroupRecord> & mReadGroups = samFileHeader->getReadGroupRecord();
    for(auto & readGroup : mReadGroups)
    {
        samplesList.emplace_back(readGroup.getReadGroupId());
    }
    NaturalLogUtils::initial();
    assert(aligner != nullptr);
    if (!modelPath.empty()) {
        mymodel.Initial(modelPath);
    }
}

Mutect2Engine::~Mutect2Engine()
{
    delete likelihoodCalculationEngine;
    delete aligner;
}


std::shared_ptr<ActivityProfileState> Mutect2Engine::isActive(AlignmentContext &context) {
	hts_pos_t pos = context.getPosition();

	std::string refName = context.getRefName();

	if (context.getReadNum() > minCallableDepth)
		callableSites++;

	char refBase = refCache->getBase(pos);

	ReadPileup tumorPileup = context.makeTumorPileup();

	std::shared_ptr<std::vector<char>> tumorAltQuals = altQuals(tumorPileup, refBase, 40);
	if (!tumorAltQuals->size())
		return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 0.0);
	double tumorLogOdds = logLikelihoodRatio(tumorPileup.size() - tumorAltQuals->size(), tumorAltQuals);

	if (tumorLogOdds < M2ArgumentCollection::getInitialLogOdds()) {
		return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 0.0);
	} else if (hasNormal() && !MTAC.genotypeGermlineSites) {
		ReadPileup normalPileup = context.makeNormalPileup();
		std::shared_ptr<std::vector<char>> normalAltQuals = altQuals(normalPileup, refBase, 40);
		int normalAltCount = normalAltQuals->size();
		double normalQualSum = 0.0;
		for (char i: *normalAltQuals) {
			normalQualSum += i;
		}
		if (normalAltCount > normalPileup.size() * 0.3 && normalQualSum > 100) {
			return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 0.0);
		}
	}
	return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 1.0);
}


std::shared_ptr<std::vector<char>> Mutect2Engine::altQuals(ReadPileup &pileup, char refBase, int pcrErrorQual) {
	std::shared_ptr<std::vector<char>> result = std::make_shared<std::vector<char>>(std::vector<char>());
	hts_pos_t pos = pileup.getPosition();

	const std::list<pileRead *> &pileupElements = pileup.getPileupElements();

	for (const pileRead *read: pileupElements) {
		PeUtils pe(read->read.get(), pos);
		int indelLength = getCurrentOrFollowingIndelLength(pe);

		if (indelLength > 0) {
			result->emplace_back(indelQual(indelLength));
		} else if (isNextToUsefulSoftClip(pe)) {
			result->emplace_back(indelQual(1));
		} else if (pe.getBase() != refBase && pe.getQual() > 6) {

			int mateStart = (!read->read->isProperlyPaired() || read->read->mateIsUnmapped()) ? INT32_MAX
			                                                                                  : read->read->getMateStart();
			bool overlapsMate = mateStart <= pos && pos < mateStart + read->read->getLength();
			result->emplace_back(
					overlapsMate ? std::min(static_cast<int>(pe.getQual()), pcrErrorQual / 2) : pe.getQual());
		}
	}
	return result;
}


char Mutect2Engine::indelQual(int indelLength) {
	return (char) std::min(30 + (indelLength - 1) * 10, 127);
}

bool Mutect2Engine::isNextToUsefulSoftClip(PeUtils &pe) {
	int pos = pe.getOffset();
	return pe.getQual() > MINIMUM_BASE_QUALITY &&
	       ((pe.isBeforeSoftClip() && pe.getBaseQuality(pos + 1) > MINIMUM_BASE_QUALITY)
	        || (pe.isAfterSoftClip() && pe.getBaseQuality(pos - 1) > MINIMUM_BASE_QUALITY));
}

int Mutect2Engine::getCurrentOrFollowingIndelLength(PeUtils &pe) {
	return pe.isDeletion() ? pe.getCurrentCigarElement()->getLength() : pe.getLengthOfImmediatelyFollowingIndel();
}

double Mutect2Engine::logLikelihoodRatio(int refCount, const std::shared_ptr<std::vector<char>> &altQuals) {
	return logLikelihoodRatio(refCount, altQuals, 1);
}

double
Mutect2Engine::logLikelihoodRatio(int nRef, const std::shared_ptr<std::vector<char>> &altQuals, int repeatFactor) {
	int nAlt = repeatFactor * altQuals->size();
	int n = nRef + nAlt;

	double fTildeRatio = std::exp(MathUtils::digamma(nRef + 1) - MathUtils::digamma(nAlt + 1));
	double betaEntropy = MathUtils::log10ToLog(
			-MathUtils::log10Factorial(n + 1) + MathUtils::log10Factorial(nAlt) + MathUtils::log10Factorial(nRef));
	double readSum = 0;
	for (char qual: *altQuals) {
		double epsilon = QualityUtils::qualToErrorProb(static_cast<uint8_t>(qual));
		double zBarAlt = (1 - epsilon) / (1 - epsilon + epsilon * fTildeRatio);
		double logEpsilon = NaturalLogUtils::qualToLogErrorProb(static_cast<uint8_t>(qual));
		double logOneMinusEpsilon = NaturalLogUtils::qualToLogProb(static_cast<uint8_t>(qual));
		readSum += zBarAlt * (logOneMinusEpsilon - logEpsilon) + MathUtils::fastBernoulliEntropy(zBarAlt);
	}
	return betaEntropy + readSum * repeatFactor;
}

bool Mutect2Engine::hasNormal() {
	return !normalSample.empty();
}

void
Mutect2Engine::fillNextAssemblyRegionWithReads(const std::shared_ptr<AssemblyRegion> &region, ReadCache &readCache) {
	std::vector<std::shared_ptr<SAMRecord>> toAdd = readCache.getReadsForRegion(*region);
	for (auto read: toAdd) {
		region->add(read);
	}
	region->sortReadsByCoordinate();
}

std::vector<std::shared_ptr<VariantContext>>
Mutect2Engine::callRegion(const std::shared_ptr<AssemblyRegion> &originalAssemblyRegion,
                          ReferenceContext &referenceContext) {
    if(originalAssemblyRegion->getStart() == 36211225)
        cout << "=============\n";

	// divide PCR qual by two in order to get the correct total qual when treating paired reads as independent
	AssemblyBasedCallerUtils::cleanOverlappingReadPairs(originalAssemblyRegion->getReads(), samplesList, normalSample, false,
	                                                    MTAC.pcrSnvQual / 2, MTAC.pcrIndelQual / 2);
	if (originalAssemblyRegion->getReads().empty())
		return {};

	removeUnmarkedDuplicates(originalAssemblyRegion);

	std::shared_ptr<AssemblyRegion> assemblyActiveRegion = AssemblyBasedCallerUtils::assemblyRegionWithWellMappedReads(
			originalAssemblyRegion, READ_QUALITY_FILTER_THRESHOLD, header);
	//if (assemblyActiveRegion->getStart() + 1 != 210416389) return {};
	std::shared_ptr<AssemblyResultSet> untrimmedAssemblyResult
			= AssemblyBasedCallerUtils::assembleReads(assemblyActiveRegion, MTAC, header, *refCache, assemblyEngine);
	std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &allVariationEvents
			= untrimmedAssemblyResult->getVariationEvents(1);
	//printVariationEvents(assemblyActiveRegion, allVariationEvents);

	std::shared_ptr<AssemblyRegionTrimmer_Result> trimmingResult
			= trimmer.trim(originalAssemblyRegion, allVariationEvents);
	if (!trimmingResult->isVariationPresent()) {
		untrimmedAssemblyResult->deleteEventMap();
		return {};
	}
	// trimmingResult->printInfo();

	std::shared_ptr<AssemblyResultSet> assemblyResult
			= trimmingResult->getNeedsTrimming() ? untrimmedAssemblyResult->trimTo(trimmingResult->getCallableRegion())
			                                     : untrimmedAssemblyResult;
	if (!assemblyResult->isisVariationPresent()) {
		untrimmedAssemblyResult->deleteEventMap();
		return {};
	}
	//assemblyResult->printSortedHaplotypes();

	std::shared_ptr<AssemblyRegion> regionForGenotyping = assemblyResult->getRegionForGenotyping();
	removeReadStubs(regionForGenotyping);

	std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> reads
			= splitReadsBySample(regionForGenotyping->getReads());
	//printReadsMap(reads);

	if (mymodel.isInitialized() && regionForGenotyping->getReads().size() > 120) {
		std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &VariationEvents
				= assemblyResult->getVariationEvents(1);
		if (!mymodel.modelRefer(reads, VariationEvents, regionForGenotyping, refCache)) {
			untrimmedAssemblyResult->deleteEventMap();
			assemblyResult->deleteEventMap();
			return {};
		}
	}

    //cerr << *originalAssemblyRegion;
    auto readLikelihoods = likelihoodCalculationEngine->computeReadLikelihoods(*assemblyResult, samplesList, *reads);
    readLikelihoods->switchToNaturalLog();

    shared_ptr<unordered_map<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>> readRealignments = AssemblyBasedCallerUtils::realignReadsToTheirBestHaplotype(*readLikelihoods, assemblyResult->getReferenceHaplotype(), assemblyResult->getPaddedReferenceLoc(), aligner);
    readLikelihoods->changeEvidence(readRealignments);

    CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(readLikelihoods, *assemblyResult, referenceContext, *regionForGenotyping->getSpan(), header);

    //---print the called variant
    std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> calls = calledHaplotypes.getCalls();
    for(auto& call : *calls)
    {
        cerr << call->getContig() << " " << call->getStart() << " " << call->getEnd() << "\n";
    }

	// Break the circular reference of pointer
	untrimmedAssemblyResult->deleteEventMap();
	assemblyResult->deleteEventMap();
	delete readLikelihoods;
    return  {allVariationEvents.begin(), allVariationEvents.end()};
}

void Mutect2Engine::removeUnmarkedDuplicates(const std::shared_ptr<AssemblyRegion> &assemblyRegion) {
	// The order of reads may causes a diffrent result. at present, it seems to have little impact
	// Therefore, we also sorted the gatk version's reads here in advance to facilitate program comparison
	std::map<std::pair<std::string, int>, std::vector<std::shared_ptr<SAMRecord>>> possibleDuplicates;
	for (auto &read: assemblyRegion->getReads()) {
		if (read->isPaired() && !read->mateIsUnmapped() &&
		    (read->getMateContig() != read->getContig() ||
		     (std::abs(read->getFragmentLength()) > HUGE_FRAGMENT_LENGTH))) {
			std::string sampleName = read->getGroup() == 0 ? "normal" : "tumor";
			std::pair<std::string, int> toAdd{sampleName, (read->isFirstOfPair() ? 1 : -1) * read->getUnclippedStart()};
			if (possibleDuplicates.find(toAdd) != possibleDuplicates.end()) {
				possibleDuplicates.find(toAdd)->second.emplace_back(read);
			} else {
				possibleDuplicates.insert({toAdd, {read}});
			}
		} else continue;
	}

	std::vector<std::shared_ptr<SAMRecord>> duplicates;

	for (const auto &group: possibleDuplicates) {
		std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> readsByContig;
		for (const auto &read: group.second) {
			if (readsByContig.find(read->getMateContig()) != readsByContig.end()) {
				readsByContig.find(read->getMateContig())->second.emplace_back(read);
			} else {
				readsByContig.insert({read->getMateContig(), {read}});
			}
		}

		for (const auto &contigReads: readsByContig) {
			int skip = contigReads.second.size() > group.second.size() / 2 ? 1 : 0;
			for (const std::shared_ptr<SAMRecord> &read: contigReads.second) {
				if (skip > 0) {
					skip--;
					continue;
				} else {
					duplicates.emplace_back(read);
				}
			}
		}
	}
	assemblyRegion->removeAll(duplicates);
}

void Mutect2Engine::removeReadStubs(const std::shared_ptr<AssemblyRegion> &assemblyRegion) {
	std::vector<std::shared_ptr<SAMRecord>> readStubs;
	std::vector<std::shared_ptr<SAMRecord>> &reads = assemblyRegion->getReads();
	readStubs.reserve(reads.size());
	for (const std::shared_ptr<SAMRecord> &read: reads) {
		if (read->getLength() < AssemblyBasedCallerUtils::MINIMUM_READ_LENGTH_AFTER_TRIMMING) {
			readStubs.emplace_back(read);
		}
	}
	assemblyRegion->removeAll(readStubs);
}

std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>
Mutect2Engine::splitReadsBySample(const std::vector<std::shared_ptr<SAMRecord>> &reads) {
	return AssemblyBasedCallerUtils::splitReadsBySample(samplesList, normalSample, reads);
}

void Mutect2Engine::setReferenceCache(ReferenceCache *cache) {
	assert(cache != nullptr);
	refCache = cache;
}

void Mutect2Engine::printVariationEvents(const std::shared_ptr<AssemblyRegion> &region,
                                         const std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &ves) {
	std::cout << "region: " << region->getStart() + 1 << " " << region->getEnd() + 1 << std::endl;
	std::cout << "allVariationEvents " << ves.size() << std::endl;
	for (const auto &ve: ves) {
		std::cout << ve->getStart() + 1 << " " << ve->getEnd() + 1 << " " << ve->getTypeString() << "\t";
		for (const auto &alt: ve->getAlternateAlleles()) {
			std::cout << ve->getReference()->getBaseString() << "==>" << alt->getBaseString() << " ";
		}
		std::cout << std::endl;
	}
}

void
Mutect2Engine::printReadsMap(const std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>& reads) {
	for (const auto &keyValuePair: *reads) {
		std::cout << keyValuePair.first << " " << keyValuePair.second.size() << std::endl;
		for (const auto &read: keyValuePair.second) {
			std::cout << read->getName() << "\t" << read->getStart() + 1 << " " << read->getEnd() + 1 << "\t";
			for (const auto &ce: read->getCigarElements()) {
				std::cout << ce.getLength() << CigarOperatorUtils::enumToCharacter(ce.getOperator());
			}
			std::cout << std::endl;
		}
	}
}
