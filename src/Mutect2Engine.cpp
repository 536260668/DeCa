//
// Created by lhh on 10/19/21.
//

#include "Mutect2Engine.h"
#include "cigar/Cigar.h"
#include "utils/PeUtils.h"


Mutect2Engine::Mutect2Engine(M2ArgumentCollection & MTAC, char * ref):minCallableDepth(MTAC.callableDepth),
                                                            normalSamples(MTAC.normalSamples) ,callableSites(0), refCache(ref)
{

}

double Mutect2Engine::logLikelihoodRatio(int refCount, std::vector<char> altQuals)
{

}

ActivityProfileState Mutect2Engine::isActive(AlignmentContext* context, ReferenceContext& ref)
{
    if(context == nullptr || context->getReadNum() == 0)
        return ActivityProfileState(context->getRefName(), context->getPosition(), 0.0);

    hts_pos_t pos = context->getPosition();

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

    altQuals(tumorPileup, refBase, 40);
    return ActivityProfileState(refName, pos, 0.0);
}

// TODO: finish this method 2021.11.1
std::vector<char> Mutect2Engine::altQuals(ReadPileup &pileup, char refBase, int pcrErrorQual)
{
    std::vector<char> result;
    hts_pos_t pos = pileup.getPosition();

    std::vector<bam1_t*> & pileupElements = pileup.getPileupElements();
    for(bam1_t* read : pileupElements)
    {
        PeUtils pe(read, pos);
        int indelLength = getCurrentOrFollowingIndelLength(read, pos);
        if(indelLength > 0) {
            result.emplace_back(indelQual(indelLength));
        } else if (isNextToUsefulSoftClip(read, pos)) {

        }
    }
}


char Mutect2Engine::indelQual(int indelLength) {
    return (char) std::min(30 + (indelLength - 1) * 10, 127);
}

bool Mutect2Engine::isNextToUsefulSoftClip(bam1_t *pe, int pos) {
    hts_pos_t offset = pe->core.pos - pos;
    uint8_t * qual = bam_get_qual(pe);

}

int Mutect2Engine::getCurrentOrFollowingIndelLength(PeUtils & pe) {
    return pe.isDeletion() ? pe.getCurrentCigarElement().getLength() : pe.getLengthOfImmediatelyFollowingIndel();
}
