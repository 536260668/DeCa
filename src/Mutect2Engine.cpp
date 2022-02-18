//
// Created by lhh on 10/19/21.
//

#include "Mutect2Engine.h"
#include "cigar/Cigar.h"
#include "utils/PeUtils.h"
#include "samtools/SAMRecord.h"
#include <cmath>
#include <utility>
#include "MathUtils.h"
#include "QualityUtils.h"
#include "NaturalLogUtils.h"
#include "AssemblyResultSet.h"
#include "haplotypecaller/AssemblyBasedCallerUtils.h"

Mutect2Engine::Mutect2Engine(M2ArgumentCollection & MTAC, char * ref, SAMFileHeader* samFileHeader):MATC(MTAC), minCallableDepth(MTAC.callableDepth),
                                                            normalSamples(MTAC.normalSamples) ,callableSites(0), refCache(ref, header), header(samFileHeader),
                                                                                                    assemblyEngine(0, 1, 128, false, false, {10, 25}),
                                                                                                    trimmer(&assemblerArgs, &header->getSequenceDictionary(), false,
                                                                                                            false)
{

}



ActivityProfileState Mutect2Engine::isActive(AlignmentContext& context, ReferenceContext& ref)
{
    hts_pos_t pos = context.getPosition();
//    if(pos == 13272)
//        std::cout << "hello" << std::endl;

    std::string refName = context.getRefName();

    if(context.getReadNum() > minCallableDepth)
        callableSites++;

    char refBase = refCache.getBase(pos);
    //std::cout << "refBase: " << refBase << std::endl;
    // TODO: divide the pileup to tumor pileup and normal pileup

    ReadPileup tumorPileup = context.makeTumorPileup();
    // TODO: calculate the activeProb

    std::vector<char> tumorAltQuals = altQuals(tumorPileup, refBase, 40);
//    if(pos == 13272) {
//        std::cout << "hello";
//    }
    double tumorLogOdds = logLikelihoodRatio(tumorPileup.size() - tumorAltQuals.size(), tumorAltQuals);

    if(tumorLogOdds < M2ArgumentCollection::getInitialLogOdds()) {
        return {refName.c_str(), pos, 0.0};
    } else if (hasNormal() && !MATC.genotypeGermlineSites) {
        ReadPileup normalPileup = context.makeNormalPileup();
        std::vector<char> normalAltQuals = altQuals(tumorPileup, refBase, 40);
        int normalAltCount = normalAltQuals.size();
        double normalQualSum = 0.0;
        for (char i : normalAltQuals) {
            normalQualSum += i;
        }
        if(normalAltCount > normalPileup.size() * 0.3 && normalQualSum > 100) {
            return {refName.c_str(), pos, 0.0};
        }
    }
    return {refName.c_str(), pos, 1.0};
}

// TODO: finish this method 2021.11.1
std::vector<char> Mutect2Engine::altQuals(ReadPileup &pileup, char refBase, int pcrErrorQual)
{
    std::vector<char> result;
    hts_pos_t pos = pileup.getPosition();

    std::vector<std::shared_ptr<SAMRecord>>  pileupElements = pileup.getPileupElements();

    for(const std::shared_ptr<SAMRecord>& read : pileupElements)
    {
//        if(pos == 13272) {
//            for (int i = 0; i < 100; i++)
//                std::cout << i << " : " << (int)bam_get_qual(read)[i] << std::endl;
//            std::cout << bam_get_qual(read) << std::endl;
//        }
        PeUtils pe(read.get(), pos);
//        uint8_t base = pe.getBase();
//        uint8_t qual = pe.getQual();
//        SAMRecord tmp(read, header);
        int indelLength = getCurrentOrFollowingIndelLength(pe);

        if(indelLength > 0) {
            result.emplace_back(indelQual(indelLength));
        } else if (isNextToUsefulSoftClip(pe, pos)) {
            result.emplace_back(indelQual(1));
        } else if (pe.getBase() != refBase && pe.getQual() > 6){

            int mateStart = (!read->isProperlyPaired() || read->mateIsUnmapped()) ? INT32_MAX : read->getMateStart();
            bool overlapsMate = mateStart <= pos && pos < mateStart + read->getLength();
            result.emplace_back(overlapsMate ? std::min(static_cast<int>(pe.getQual()), pcrErrorQual/2) : pe.getQual());
        }
    }
    return result;
}


char Mutect2Engine::indelQual(int indelLength) {
    return (char) std::min(30 + (indelLength - 1) * 10, 127);
}

bool Mutect2Engine::isNextToUsefulSoftClip(PeUtils & pe, int pos) {
    return pe.getQual() > MINIMUM_BASE_QUALITY &&
            ((pe.isBeforeSoftClip() && pe.getBaseQuality(pos + 1) > MINIMUM_BASE_QUALITY)
            || (pe.isAfterSoftClip() && pe.getBaseQuality(pos - 1) > MINIMUM_BASE_QUALITY));
}

int Mutect2Engine::getCurrentOrFollowingIndelLength(PeUtils & pe) {
    return pe.isDeletion() ? pe.getCurrentCigarElement().getLength() : pe.getLengthOfImmediatelyFollowingIndel();
}

double Mutect2Engine::logLikelihoodRatio(int refCount, std::vector<char> &altQuals) {
    return logLikelihoodRatio(refCount, altQuals, 1);
}

double Mutect2Engine::logLikelihoodRatio(int nRef, std::vector<char> &altQuals, int repeatFactor) {
    int nAlt = repeatFactor * altQuals.size();
    int n = nRef + nAlt;

    double fTildeRatio = std::exp(MathUtils::digamma(nRef + 1) - MathUtils::digamma(nAlt + 1));
    double betaEntropy = MathUtils::log10ToLog(-MathUtils::log10Factorial(n+1) + MathUtils::log10Factorial(nAlt) + MathUtils::log10Factorial(nRef));
    double readSum = 0;
    for(char qual : altQuals) {
        double epsilon = QualityUtils::qualToErrorProb(static_cast<uint8_t>(qual));
        double zBarAlt = (1 - epsilon) / (1 - epsilon + epsilon * fTildeRatio);
        double logEpsilon = NaturalLogUtils::qualToLogErrorProb(static_cast<uint8_t>(qual));
        double logOneMinusEpsilon = NaturalLogUtils::qualToLogProb(static_cast<uint8_t>(qual));
        readSum += zBarAlt * (logOneMinusEpsilon - logEpsilon) + MathUtils::fastBernoulliEntropy(zBarAlt);
    }
    return betaEntropy + readSum * repeatFactor;
}

bool Mutect2Engine::hasNormal() {
    return !normalSamples.empty();
}

void Mutect2Engine::fillNextAssemblyRegionWithReads(const std::shared_ptr<AssemblyRegion>& region, ReadCache &readCache) {
    std::vector<std::shared_ptr<SAMRecord>> toAdd = readCache.getReadsForRegion(*region);
    region->setRead(toAdd);
}

std::vector<std::shared_ptr<VariantContext>>
Mutect2Engine::callRegion(std::shared_ptr<AssemblyRegion> originalAssemblyRegion, ReferenceContext &referenceContext) {
    if(originalAssemblyRegion->getStart() == 1017765) {
        for(const std::shared_ptr<SAMRecord>& read : originalAssemblyRegion->getReads()) {
            std::cout << read->getName() << " : " << read->getStart() + 1 << "~" << read->getEnd() + 1 << std::endl;
        }
    }
    removeUnmarkedDuplicates(originalAssemblyRegion);
    if(originalAssemblyRegion->getReads().size() == 0)
        return {};
    std::shared_ptr<AssemblyResultSet> untrimmedAssemblyResult = AssemblyBasedCallerUtils::assembleReads(std::move(originalAssemblyRegion), MATC, header, refCache, assemblyEngine);
    return {};
    //std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & allVariationEvents = untrimmedAssemblyResult->getVariationEvents(1);
    //std::shared_ptr<AssemblyRegionTrimmer_Result> trimmingResult = trimmer.trim(originalAssemblyRegion, allVariationEvents);
    //return  {allVariationEvents.begin(), allVariationEvents.end()};
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
