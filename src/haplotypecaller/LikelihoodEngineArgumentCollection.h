//
// Created by lhh on 4/21/22.
//

#ifndef MUTECT2CPP_MASTER_LIKELIHOODENGINEARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_LIKELIHOODENGINEARGUMENTCOLLECTION_H

#include "utils/pairhmm/PairHMM.h"
#include "PairHMMLikelihoodCalculationEngine.h"
#include "PairHMMNativeArgumentCollection.h"

class LikelihoodEngineArgumentCollection {
public:
    /**
     * Bases with a quality below this threshold will reduced to the minimum usable qualiy score (6).
     */
    char BASE_QUALITY_SCORE_THRESHOLD = PairHMM::BASE_QUALITY_SCORE_THRESHOLD;

    int gcpHMM = 10;

    /**
     * When calculating the likelihood of variants, we can try to correct for PCR errors that cause indel artifacts.
     * The correction is based on the reference context, and acts specifically around repetitive sequences that tend
     * to cause PCR errors). The variant likelihoods are penalized in increasing scale as the context around a
     * putative indel is more repetitive (e.g. long homopolymer). The correction can be disabling by specifying
     * '-pcrModel NONE'; in that case the default base insertion/deletion qualities will be used (or taken from the
     * read if generated through the BaseRecalibrator). <b>VERY IMPORTANT: when using PCR-free sequencing data we
     * definitely recommend setting this argument to NONE</b>.
     */
    PCRErrorModel pcrErrorModel = CONSERVATIVE;

    /**
     * The phredScaledGlobalReadMismappingRate reflects the average global mismapping rate of all reads, regardless of their
     * mapping quality.  This term effects the probability that a read originated from the reference haplotype, regardless of
     * its edit distance from the reference, in that the read could have originated from the reference haplotype but
     * from another location in the genome.  Suppose a read has many mismatches from the reference, say like 5, but
     * has a very high mapping quality of 60.  Without this parameter, the read would contribute 5 * Q30 evidence
     * in favor of its 5 mismatch haplotype compared to reference, potentially enough to make a call off that single
     * read for all of these events.  With this parameter set to Q30, though, the maximum evidence against any haplotype
     * that this (and any) read could contribute is Q30.
     *
     * Set this term to any negative number to turn off the global mapping rate.
     */
    int phredScaledGlobalReadMismappingRate = 45;

    PairHMMNativeArgumentCollection pairHMMNativeArgs;
};


#endif //MUTECT2CPP_MASTER_LIKELIHOODENGINEARGUMENTCOLLECTION_H
