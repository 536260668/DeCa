//
// Created by 梦想家xixi on 2021/12/8.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H
#define MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H

#include "haplotype/Haplotype.h"
#include "AssemblyRegion.h"
#include "ReferenceCache.h"
#include "Mutect2/AssemblyResultSet.h"
#include "M2ArgumentCollection.h"
#include "ReadThreadingAssembler.h"
#include "PairHMMLikelihoodCalculationEngine.h"
#include "LikelihoodEngineArgumentCollection.h"

static const std::string NORMAL = "normal";
static const std::string TUMOR = "tumor";

class AssemblyBasedCallerUtils {
public:
    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    static std::shared_ptr<Haplotype> createReferenceHaplotype(const std::shared_ptr<AssemblyRegion> & region, const std::shared_ptr<SimpleInterval> &referencePadding, ReferenceCache & cache);

    static std::shared_ptr<AssemblyResultSet> assembleReads(const std::shared_ptr<AssemblyRegion>& region, M2ArgumentCollection & argumentCollection, SAMFileHeader* header, ReferenceCache & cache, ReadThreadingAssembler& assemblyEngine);

    static const int REFERENCE_PADDING_FOR_ASSEMBLY = 500;

    static const int MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;

    static std::shared_ptr<SimpleInterval> getPaddedReferenceLoc(const std::shared_ptr<AssemblyRegion>& region, int referencePadding, SAMFileHeader* header);

    static void finalizeRegion(const std::shared_ptr<AssemblyRegion>& region, bool errorCorrectReads, bool dontUseSoftClippedBases, uint8_t minTailQuality, SAMFileHeader* header, bool correctOverlappingBaseQualities);

    static std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> splitReadsBySample(const std::vector<std::shared_ptr<SAMRecord>> & reads);

    /**
     * Instantiates the appropriate likelihood calculation engine.
     */
    static PairHMMLikelihoodCalculationEngine * createLikelihoodCalculationEngine(LikelihoodEngineArgumentCollection& likelihoodArgs);

    /**
     *  Modify base qualities when paired reads overlap to account for the possibility of PCR error.
     *
     *  Overlapping mates provded independent evidence as far as sequencing error is concerned, but their PCR errors
     *  are correlated.  The base qualities are thus limited by the sequencing base quality as well as half of the PCR
     *  quality.  We use half of the PCR quality because downstream we treat read pairs as independent, and summing two halves
     *  effectively gives the PCR quality of the pairs when taken together.
     *
     * @param reads the list of reads to consider
     * @param samplesList   list of samples|
     * @param readsHeader   bam header of reads' source
     * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
     * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR
     * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR
     */
    static void cleanOverlappingReadPairs(vector<shared_ptr<SAMRecord>>& reads, const string& sample, bool setConflictingToZero, int halfOfPcrSnvQual = 0, int halfOfPcrIndelQual = 0);

    // create the assembly using just high quality reads (eg Q20 or higher).  We may want to use lower
    // quality reads in the PairHMM downstream, so we can't use a ReadFilter
    static std::shared_ptr<AssemblyRegion> assemblyRegionWithWellMappedReads(const std::shared_ptr<AssemblyRegion>& originalAssemblyRegion, int minMappingQuality, SAMFileHeader * header);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H
