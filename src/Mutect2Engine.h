//
// The core algorithm of Mutect2 is here
// Created by lhh on 10/19/21.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2ENGINE_H
#define MUTECT2CPP_MASTER_MUTECT2ENGINE_H

#include <vector>
#include "ActivityProfileState.h"
#include "M2ArgumentCollection.h"
#include "ReferenceCache.h"
#include "engine/ReferenceContext.h"
#include "utils/PeUtils.h"
#include "samtools/SAMFileHeader.h"
#include "engine/AlignmentContext.h"
#include "utils/ReadPileup.h"
#include "ReadCache.h"

class Mutect2Engine {
private:
    int minCallableDepth;
    std::set<std::string> & normalSamples;
    SAMFileHeader * header;
    M2ArgumentCollection & MATC;



    std::vector<char> altQuals(ReadPileup & pileup, char refBase, int pcrErrorQual);

    static int getCurrentOrFollowingIndelLength(PeUtils & pe);

    static char indelQual(int indelLength);

    static bool isNextToUsefulSoftClip(PeUtils & pe, int pos);

    /**
     * this implement the isActive() algorithm described in docs/mutect/mutect.pdf
     * the multiplicative factor is for the special case where we pass a singleton list
     * of alt quals and want to duplicate that alt qual over multiple reads
     */
    static double logLikelihoodRatio(int refCount, std::vector<char> & altQuals);

    static double logLikelihoodRatio(int nRef, std::vector<char> & altQuals, int repeatFactor);

    bool hasNormal();
public:
    const static int READ_QUALITY_FILTER_THRESHOLD = 20;
    const static int MIN_READ_LENGTH = 30;
    const static int MINIMUM_BASE_QUALITY = 6;

    int callableSites;  // in GATK4, this variable is a MutableInt class object
    ReferenceCache refCache;

    Mutect2Engine(M2ArgumentCollection & MTAC, char* ref, SAMFileHeader*);

    ActivityProfileState isActive(AlignmentContext& context, ReferenceContext & referenceContext);

    static void fillNextAssemblyRegionWithReads(AssemblyRegion & region, ReadCache & readCache);
};


#endif //MUTECT2CPP_MASTER_MUTECT2ENGINE_H
