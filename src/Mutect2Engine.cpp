//
// Created by lhh on 10/19/21.
//

#include "Mutect2Engine.h"
#include "cigar/Cigar.h"
#include "utils/PeUtils.h"
#include "samtools/SAMRecord.h"
#include <cmath>
#include "MathUtils.h"
#include "QualityUtils.h"
#include "NaturalLogUtils.h"

Mutect2Engine::Mutect2Engine(M2ArgumentCollection & MTAC, char * ref, SAMFileHeader* samFileHeader):MATC(MTAC), minCallableDepth(MTAC.callableDepth),
                                                            normalSamples(MTAC.normalSamples) ,callableSites(0), refCache(ref), header(samFileHeader)
{

}



ActivityProfileState Mutect2Engine::isActive(AlignmentContext* context, ReferenceContext& ref)
{
    if(context == nullptr || context->getReadNum() == 0)
        return ActivityProfileState(context->getRefName(), context->getPosition(), 0.0);

    hts_pos_t pos = context->getPosition();

//    if(pos == 13272)
//        std::cout << "hello" << std::endl;

    const char * refName = context->getRefName();

    if(context->getReadNum() > minCallableDepth)
        callableSites++;



    // TODO: add ApplyBQSR for reads


    if (strcmp(refName, refCache.getContig()) != 0 )
    {
        refCache.clear();
        /*if(DEBUG_OPEN)
        {
            cout << refName << " " << refCache.getContig() << endl;
        }*/
        refCache.fill(refName, pos);
    }
    char refBase = refCache.getBase(pos);
    //std::cout << "refBase: " << refBase << std::endl;
    // TODO: divide the pileup to tumor pileup and normal pileup

    ReadPileup tumorPileup = context->makeTumorPileup(normalSamples);
    // TODO: calculate the activeProb

    std::vector<char> tumorAltQuals = altQuals(tumorPileup, refBase, 40);
    double tumorLogOdds = logLikelihoodRatio(tumorPileup.size() - tumorAltQuals.size(), tumorAltQuals);
    if(tumorLogOdds < M2ArgumentCollection::getInitialLogOdds()) {
        return {refName, pos, 0.0};
    } else if (hasNormal() && !MATC.genotypeGermlineSites) {
        ReadPileup normalPileup = context->makeNormalPileup(normalSamples);
        std::vector<char> normalAltQuals = altQuals(tumorPileup, refBase, 40);
        int normalAltCount = normalAltQuals.size();
        double normalQualSum = 0.0;
        for (char i : normalAltQuals) {
            normalQualSum += i;
        }
        if(normalAltCount > normalPileup.size() * 0.3 && normalQualSum > 100) {
            return {refName, pos, 0.0};
        }
    }

    return ActivityProfileState(refName, pos, 0.0);
}

// TODO: finish this method 2021.11.1
std::vector<char> Mutect2Engine::altQuals(ReadPileup &pileup, char refBase, int pcrErrorQual)
{
    std::vector<char> result;
    hts_pos_t pos = pileup.getPosition();

    std::vector<bam1_t*>  pileupElements = pileup.getPileupElements();

    for(bam1_t* read : pileupElements)
    {
//        if(pos == 13272) {
//            for (int i = 0; i < 100; i++)
//                std::cout << i << " : " << (int)bam_get_qual(read)[i] << std::endl;
//            std::cout << bam_get_qual(read) << std::endl;
//        }
        PeUtils pe(read, pos);
//        uint8_t base = pe.getBase();
//        uint8_t qual = pe.getQual();
//        SAMRecord tmp(read, header);
        int indelLength = getCurrentOrFollowingIndelLength(pe);

        if(indelLength > 0) {
            result.emplace_back(indelQual(indelLength));
        } else if (isNextToUsefulSoftClip(pe, pos)) {
            result.emplace_back(indelQual(1));
        } else if (pe.getBase() != refBase && pe.getQual() > 6){

            SAMRecord samRecord(read, header);
            int mateStart = (!samRecord.isProperlyPaired() || samRecord.mateIsUnmapped()) ? INT32_MAX : samRecord.getMateStart();
            bool overlapsMate = mateStart <= pos && pos < mateStart + samRecord.getLength();
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
