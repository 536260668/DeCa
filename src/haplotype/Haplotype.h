//
// Created by 梦想家xixi on 2021/11/9.
//

#ifndef MUTECT2CPP_MASTER_HAPLOTYPE_H
#define MUTECT2CPP_MASTER_HAPLOTYPE_H


#include "Allele.h"
#include "Locatable.h"
#include "cigar/Cigar.h"

class Haplotype : public Allele{
private:
    Locatable* genomeLocation;
    Cigar cigar;
    int alignmentStartHapwrtRef;
    double score;
    static uint8_t * copyArray(uint8_t * base, int length);

public:
    /**
     * Main constructor
     *
     * @param bases a non-null array of bases
     * @param isRef is this the reference haplotype?
     */
    Haplotype(uint8_t *bases, int length, bool isRef);

    /**
     * Create a new non-ref haplotype
     *
     * @param bases a non-null array of bases
     */
     Haplotype(uint8_t * bases, int length);

    /**
    * Create a new haplotype with bases
    *
    * Requires bases.length == cigar.getReadLength()
    *
    * @param bases a non-null array of bases
    * @param isRef is this the reference haplotype?
    * @param alignmentStartHapwrtRef offset of this haplotype w.r.t. the reference
    * @param cigar the cigar that maps this haplotype to the reference sequence
    */
    Haplotype(uint8_t* bases, bool isRef, int length, int alignmentStartHapwrtRef, Cigar * cigar);

    /**
    * Set the cigar of this haplotype to cigar.
    *
    * Note that this function consolidates the cigar, so that 1M1M1I1M1M => 2M1I2M
    *
    * @param cigar a cigar whose readLength == length()
    */
     void setCigar(Cigar *cigar);
};


#endif //MUTECT2CPP_MASTER_HAPLOTYPE_H
