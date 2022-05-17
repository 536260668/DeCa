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

Mutect2Engine::Mutect2Engine(M2ArgumentCollection & MTAC, char * ref, SAMFileHeader* samFileHeader):MTAC(MTAC), minCallableDepth(MTAC.callableDepth),
                                                            normalSample(MTAC.normalSample) ,callableSites(0), refCache(nullptr) ,header(samFileHeader),
                                                                                                    assemblyEngine(0, 1, 128, false, false, {10, 25}),
                                                                                                    likelihoodCalculationEngine(AssemblyBasedCallerUtils::createLikelihoodCalculationEngine(MTAC.likelihoodArgs)),
                                                                                                    trimmer(&assemblerArgs, &header->getSequenceDictionary(), false,
                                                                                                            false)
{
    std::vector<SAMReadGroupRecord> & mReadGroups = samFileHeader->getReadGroupRecord();
    for(auto & readGroup : mReadGroups)
    {
        samplesList.emplace_back(readGroup.getReadGroupId());
    }
    NaturalLogUtils::initial();
}

Mutect2Engine::~Mutect2Engine()
{
    delete likelihoodCalculationEngine;
}


std::shared_ptr<ActivityProfileState> Mutect2Engine::isActive(AlignmentContext& context)
{
    hts_pos_t pos = context.getPosition();

    std::string refName = context.getRefName();

    if(context.getReadNum() > minCallableDepth)
        callableSites++;

    char refBase = refCache->getBase(pos);

    ReadPileup tumorPileup = context.makeTumorPileup();

    std::shared_ptr<std::vector<char>> tumorAltQuals = altQuals(tumorPileup, refBase, 40);
    if(!tumorAltQuals->size())
        return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 0.0);
    double tumorLogOdds = logLikelihoodRatio(tumorPileup.size() - tumorAltQuals->size(), tumorAltQuals);

    if(tumorLogOdds < M2ArgumentCollection::getInitialLogOdds()) {
        return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 0.0);
    } else if (hasNormal() && !MTAC.genotypeGermlineSites) {
        ReadPileup normalPileup = context.makeNormalPileup();
        std::shared_ptr<std::vector<char>> normalAltQuals = altQuals(normalPileup, refBase, 40);
        int normalAltCount = normalAltQuals->size();
        double normalQualSum = 0.0;
        for (char i : *normalAltQuals) {
            normalQualSum += i;
        }
        if(normalAltCount > normalPileup.size() * 0.3 && normalQualSum > 100) {
            return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 0.0);
        }
    }
    return std::make_shared<ActivityProfileState>(refName.c_str(), pos, 1.0);
}


std::shared_ptr<std::vector<char>> Mutect2Engine::altQuals(ReadPileup &pileup, char refBase, int pcrErrorQual)
{
    std::shared_ptr<std::vector<char>> result = std::make_shared<std::vector<char>>(std::vector<char>());
    hts_pos_t pos = pileup.getPosition();

    const std::list<pileRead*>& pileupElements = pileup.getPileupElements();

    for(const pileRead* read : pileupElements)
    {
        PeUtils pe(read->read.get(), pos);
        int indelLength = getCurrentOrFollowingIndelLength(pe);

        if(indelLength > 0) {
            result->emplace_back(indelQual(indelLength));
        } else if (isNextToUsefulSoftClip(pe)) {
            result->emplace_back(indelQual(1));
        } else if (pe.getBase() != refBase && pe.getQual() > 6){

            int mateStart = (!read->read->isProperlyPaired() || read->read->mateIsUnmapped()) ? INT32_MAX : read->read->getMateStart();
            bool overlapsMate = mateStart <= pos && pos < mateStart + read->read->getLength();
            result->emplace_back(overlapsMate ? std::min(static_cast<int>(pe.getQual()), pcrErrorQual/2) : pe.getQual());
        }
    }
    return result;
}


char Mutect2Engine::indelQual(int indelLength) {
    return (char) std::min(30 + (indelLength - 1) * 10, 127);
}

bool Mutect2Engine::isNextToUsefulSoftClip(PeUtils & pe) {
    int pos = pe.getOffset();
    return pe.getQual() > MINIMUM_BASE_QUALITY &&
            ((pe.isBeforeSoftClip() && pe.getBaseQuality(pos + 1) > MINIMUM_BASE_QUALITY)
            || (pe.isAfterSoftClip() && pe.getBaseQuality(pos - 1) > MINIMUM_BASE_QUALITY));
}

int Mutect2Engine::getCurrentOrFollowingIndelLength(PeUtils & pe) {
    return pe.isDeletion() ? pe.getCurrentCigarElement()->getLength() : pe.getLengthOfImmediatelyFollowingIndel();
}

double Mutect2Engine::logLikelihoodRatio(int refCount, const std::shared_ptr<std::vector<char>> &altQuals) {
    return logLikelihoodRatio(refCount, altQuals, 1);
}

double Mutect2Engine::logLikelihoodRatio(int nRef, const std::shared_ptr<std::vector<char>> &altQuals, int repeatFactor) {
    int nAlt = repeatFactor * altQuals->size();
    int n = nRef + nAlt;

    double fTildeRatio = std::exp(MathUtils::digamma(nRef + 1) - MathUtils::digamma(nAlt + 1));
    double betaEntropy = MathUtils::log10ToLog(-MathUtils::log10Factorial(n+1) + MathUtils::log10Factorial(nAlt) + MathUtils::log10Factorial(nRef));
    double readSum = 0;
    for(char qual : *altQuals) {
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

void Mutect2Engine::fillNextAssemblyRegionWithReads(const std::shared_ptr<AssemblyRegion>& region, ReadCache &readCache) {
    std::vector<std::shared_ptr<SAMRecord>> toAdd = readCache.getReadsForRegion(*region);
	for (auto read : toAdd){
		region->add(read);
	}
	region->sortReadsByCoordinate();
}

std::vector<std::shared_ptr<VariantContext>>
Mutect2Engine::callRegion(const std::shared_ptr<AssemblyRegion>& originalAssemblyRegion, ReferenceContext &referenceContext) {
//    if(originalAssemblyRegion->getStart() == 359408) {
////        for(const std::shared_ptr<SAMRecord>& read : originalAssemblyRegion->getReads()) {
////            std::cout << read->getName() << " : " << read->getStart() + 1 << "~" << read->getEnd() + 1 << std::endl;
////        }
//        std::cout << "hello" << std::endl;
//    }

    // divide PCR qual by two in order to get the correct total qual when treating paired reads as independent
    AssemblyBasedCallerUtils::cleanOverlappingReadPairs(originalAssemblyRegion->getReads(), normalSample, false, MTAC.pcrSnvQual/2, MTAC.pcrIndelQual/2);

    removeUnmarkedDuplicates(originalAssemblyRegion);
    if(originalAssemblyRegion->getReads().empty())
        return {};

    auto assemblyActiveRegion = AssemblyBasedCallerUtils::assemblyRegionWithWellMappedReads(originalAssemblyRegion, READ_QUALITY_FILTER_THRESHOLD, header);
    std::shared_ptr<AssemblyResultSet> untrimmedAssemblyResult = AssemblyBasedCallerUtils::assembleReads(assemblyActiveRegion, MTAC, header, *refCache, assemblyEngine);
    std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & allVariationEvents = untrimmedAssemblyResult->getVariationEvents(1);

    std::shared_ptr<AssemblyRegionTrimmer_Result> trimmingResult = trimmer.trim(originalAssemblyRegion, allVariationEvents);
    if(!trimmingResult->isVariationPresent()) {
	    untrimmedAssemblyResult->deleteEventMap();
        return {};
    }

    std::shared_ptr<AssemblyResultSet> assemblyResult = trimmingResult->getNeedsTrimming() ? untrimmedAssemblyResult->trimTo(trimmingResult->getCallableRegion()) : untrimmedAssemblyResult;
    if(!assemblyResult->isisVariationPresent()) {
		untrimmedAssemblyResult->deleteEventMap();
        return {};
    }
    std::shared_ptr<AssemblyRegion> regionForGenotyping = assemblyResult->getRegionForGenotyping();
    removeReadStubs(regionForGenotyping);

    auto reads = splitReadsBySample(regionForGenotyping->getReads());

    //cerr << *originalAssemblyRegion;
    likelihoodCalculationEngine->computeReadLikelihoods(*assemblyResult, samplesList, *reads);

	// Break the circular reference of pointer
	untrimmedAssemblyResult->deleteEventMap();
	assemblyResult->deleteEventMap();
    return  {allVariationEvents.begin(), allVariationEvents.end()};
}

void Mutect2Engine::removeUnmarkedDuplicates(const std::shared_ptr<AssemblyRegion>& assemblyRegion) {
    std::map<std::pair<std::string, int> ,std::vector<std::shared_ptr<SAMRecord>>> possibleDuplicates;
    for(const std::shared_ptr<SAMRecord>& read : assemblyRegion->getReads()) {
        if(read->isPaired() && !read->mateIsUnmapped() &&
                (read->getMateContig() != read->getContig() || (std::abs(read->getFragmentLength()) > HUGE_FRAGMENT_LENGTH))) {
            std::string sampleName = read->getGroup() == 0 ? "normal" : "tumor";
            std::pair<std::string, int> toAdd{sampleName, (read->isFirstOfPair() ? 1 : -1) * read->getUnclippedStart()};
            if(possibleDuplicates.find(toAdd) != possibleDuplicates.end()) {
                possibleDuplicates.find(toAdd)->second.emplace_back(read);
            } else {
                possibleDuplicates.insert({toAdd, {read}});
            }
        } else {
            continue;
        }
    }

    std::vector<std::shared_ptr<SAMRecord>> duplicates;

    for(std::pair<std::pair<std::string, int> ,std::vector<std::shared_ptr<SAMRecord>>> group : possibleDuplicates) {
        std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> readsByContig;
        for(const std::shared_ptr<SAMRecord>& read : group.second) {
            if(readsByContig.find(read->getMateContig()) != readsByContig.end()) {
                readsByContig.find(read->getMateContig())->second.emplace_back(read);
            } else {
                readsByContig.insert({read->getMateContig(), {read}});
            }
        }

        for(std::pair<std::string, std::vector<std::shared_ptr<SAMRecord>>> contigReads : readsByContig) {
            int skip = contigReads.second.size() > group.second.size() / 2 ? 1 : 0;
            for(const std::shared_ptr<SAMRecord>& read : contigReads.second) {
                if(skip > 0) {
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

void Mutect2Engine::removeReadStubs(const std::shared_ptr<AssemblyRegion> & assemblyRegion) {
    std::vector<std::shared_ptr<SAMRecord>> readStubs;
    std::vector<std::shared_ptr<SAMRecord>> & reads = assemblyRegion->getReads();
    readStubs.reserve(reads.size());
    for(const std::shared_ptr<SAMRecord> & read : reads) {
        if(read->getLength() < AssemblyBasedCallerUtils::MINIMUM_READ_LENGTH_AFTER_TRIMMING) {
            readStubs.emplace_back(read);
        }
    }
    assemblyRegion->removeAll(readStubs);
}

std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>
Mutect2Engine::splitReadsBySample(const std::vector<std::shared_ptr<SAMRecord>> & reads) {
    return AssemblyBasedCallerUtils::splitReadsBySample(reads);
}

void Mutect2Engine::setReferenceCache(ReferenceCache *cache)
{
    assert(cache != nullptr);
    refCache = cache;
}
