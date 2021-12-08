//
// Created by 梦想家xixi on 2021/11/9.
//

#ifndef MUTECT2CPP_MASTER_HAPLOTYPE_H
#define MUTECT2CPP_MASTER_HAPLOTYPE_H


#include "Allele.h"
#include "Locatable.h"
#include "cigar/Cigar.h"
#include "EventMap.h"

class EventMap;

class Haplotype : public Allele{
private:
    Locatable* genomeLocation;
    Cigar* cigar;
    EventMap * eventMap;
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

    Haplotype(uint8_t* bases, int length, Locatable* loc);

    /**
    * Set the cigar of this haplotype to cigar.
    *
    * Note that this function consolidates the cigar, so that 1M1M1I1M1M => 2M1I2M
    *
    * @param cigar a cigar whose readLength == length()
    */
     void setCigar(Cigar *cigar);

    /**
    * Create a new Haplotype derived from this one that exactly spans the provided location
    *
    * Note that this haplotype must have a contain a genome loc for this operation to be successful.  If no
    * GenomeLoc is contained than @throws an IllegalStateException
    *
    * Also loc must be fully contained within this Haplotype's genomeLoc.  If not an IllegalArgumentException is
    * thrown.
    *
    * @param loc a location completely contained within this Haplotype's location
    * @return a new Haplotype within only the bases spanning the provided location, or null for some reason the haplotype would be malformed if
    */
     Haplotype* trim(Locatable* loc);

    /**
    * Get the cigar for this haplotype.  Note that the cigar is guaranteed to be consolidated
    * in that multiple adjacent equal operates will have been merged
    * @return the cigar of this haplotype
    */
     Cigar* getCigar();

     void setGenomeLocation(Locatable* genomeLocation);

     void setScore(double score);

     void setAlignmentStartHapwrtRef(int alignmentStartHapwrtRef) {
         this->alignmentStartHapwrtRef = alignmentStartHapwrtRef;
     }

    int getAlignmentStartHapwrtRef() const;

     Locatable* getGenomeLocation() {return genomeLocation;}

     EventMap* getEventMap();

     void setEventMap(EventMap* eventMap);
};


#endif //MUTECT2CPP_MASTER_HAPLOTYPE_H
