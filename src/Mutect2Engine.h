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
#include "variantcontext/VariantContext.h"
#include "ReadThreadingAssembler.h"
#include "haplotypecaller/MutectReadThreadingAssemblerArgumentCollection.h"
#include "haplotypecaller/AssemblyRegionTrimmer.h"

class Mutect2Engine {
private:
    int minCallableDepth;
    std::set<std::string> & normalSamples;
    SAMFileHeader * header;
    M2ArgumentCollection & MATC;
    ReadThreadingAssembler assemblyEngine;
    MutectReadThreadingAssemblerArgumentCollection assemblerArgs;
    AssemblyRegionTrimmer trimmer;

    std::shared_ptr<std::vector<char>> altQuals(ReadPileup & pileup, char refBase, int pcrErrorQual);

    static int getCurrentOrFollowingIndelLength(PeUtils & pe);

    static char indelQual(int indelLength);

    static bool isNextToUsefulSoftClip(PeUtils & pe);

    /**
     * this implement the isActive() algorithm described in docs/mutect/mutect.pdf
     * the multiplicative factor is for the special case where we pass a singleton list
     * of alt quals and want to duplicate that alt qual over multiple reads
     */
    static double logLikelihoodRatio(int refCount, const std::shared_ptr<std::vector<char>> & altQuals);

    static double logLikelihoodRatio(int nRef, const std::shared_ptr<std::vector<char>> & altQuals, int repeatFactor);

    bool hasNormal();

    void removeUnmarkedDuplicates(const std::shared_ptr<AssemblyRegion>& assemblyRegion);

public:
    const static int READ_QUALITY_FILTER_THRESHOLD = 20;
    const static int MIN_READ_LENGTH = 30;
    const static int MINIMUM_BASE_QUALITY = 6;

    const static int HUGE_FRAGMENT_LENGTH = 1000000;

    int callableSites;  // in GATK4, this variable is a MutableInt class object
    ReferenceCache refCache;

    Mutect2Engine(M2ArgumentCollection & MTAC, char* ref, SAMFileHeader*);

    std::shared_ptr<ActivityProfileState> isActive(AlignmentContext& context, ReferenceContext & referenceContext);

    static void fillNextAssemblyRegionWithReads(const std::shared_ptr<AssemblyRegion>& region, ReadCache & readCache);

    std::vector<std::shared_ptr<VariantContext>> callRegion(std::shared_ptr<AssemblyRegion> originalAssemblyRegion, ReferenceContext & referenceContext);

};


#endif //MUTECT2CPP_MASTER_MUTECT2ENGINE_H
