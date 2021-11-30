//
// Created by 梦想家xixi on 2021/11/10.
//

#ifndef MUTECT2CPP_MASTER_ALIGNMENTUTILS_H
#define MUTECT2CPP_MASTER_ALIGNMENTUTILS_H


#include "cigar/Cigar.h"

class AlignmentUtils {
public:
    static Cigar* consolidateCigar(Cigar* c);
    static bool needsConsolidation(Cigar* c);

    /**
     * Get the byte[] from bases that cover the reference interval refStart -> refEnd given the
     * alignment of bases to the reference (basesToRefCigar) and the start offset of the bases on the reference
     *
     * refStart and refEnd are 0 based offsets that we want to obtain.  In the client code, if the reference
     * bases start at position X and you want Y -> Z, refStart should be Y - X and refEnd should be Z - X.
     *
     * If refStart or refEnd would start or end the new bases within a deletion, this function will return null
     *
     * @param bases
     * @param refStart
     * @param refEnd
     * @param basesStartOnRef where does the bases array start w.r.t. the reference start?  For example, bases[0] of
     *                        could be at refStart == 0 if basesStartOnRef == 0, but it could just as easily be at
     *                        10 (meaning bases doesn't fully span the reference), which would be indicated by basesStartOnRef == 10.
     *                        It's not trivial to eliminate this parameter because it's tied up with the cigar
     * @param basesToRefCigar the cigar that maps the bases to the reference genome
     * @return a byte[] containing the bases covering this interval, or null if we would start or end within a deletion
     */
    static uint8_t * getBasesCoveringRefInterval(int refStart, int refEnd, uint8_t* bases, int length, int basesStartOnRef, Cigar* basesToRefCigar);

    static Cigar * trimCigarByReference(Cigar * cigar, int start, int end);

    /**
     * Does cigar start or end with a deletion operation?
     *
     * @param cigar a non-null cigar to test
     * @return true if the first or last operator of cigar is a D
     */
    static bool startsOrEndsWithInsertionOrDeletion(Cigar * cigar);

    /**
     * Removing a trailing deletion from the incoming cigar if present
     *
     * @param c the cigar we want to update
     * @return a non-null Cigar
     */
    static Cigar* removeTrailingDeletions(Cigar* c);

    static Cigar* trimCigarByBases(Cigar* cigar, int start, int end);

    static Cigar* leftAlignSingleIndel(Cigar* cigar, uint8_t* refSeq, int refLength, uint8_t* readSeq, int readLength, int refIndex, int readIndex, bool cleanupCigar);

    static Cigar* leftAlignSingleIndel(Cigar* cigar, uint8_t* refSeq, int refLength, uint8_t* readSeq, int readLength, int refIndex, int readIndex, int leftmostAllowedAlignment, bool cleanupCigar1);

    static Cigar* cleanUpCigar(Cigar* c);

private:
    static Cigar* trimCigar(Cigar * cigar, int start, int end, bool byReference);

    static void ensureLeftAlignmentHasGoodArguments(Cigar* cigar, uint8_t* refSeq, uint8_t* readSeq, int refIndex, int readIndex);

    static uint8_t * createIndelString(Cigar* cigar, int indexOfIndel, uint8_t * refSeq, int refLength, uint8_t * readSeq, int readLength, int refIndex, int readIndex, int &);

    static Cigar* moveCigarLeft(Cigar* cigar, int indexOfIndel);

protected:
    /**
     * Helper function for trimCigar that adds cigar elements (of total length X) of elt.op to dest for
     * X bases that fall between start and end, where the last position of the base is pos.
     *
     * The primary use of this function is to create a new cigar element list that contains only
     * elements that occur between start and end bases in an initial cigar.
     *
     * Note that this function may return multiple cigar elements (1M1M etc) that are best consolidated
     * after the fact into a single simpler representation.
     *
     * @param dest we will append our cigar elements to this list
     * @param pos the position (0 indexed) where elt started
     * @param start only include bases that occur >= this position
     * @param end only include bases that occur <= this position
     * @param elt the element we are slicing down
     * @return the position after we've traversed all elt.length bases of elt
     */
    static int addCigarElements(std::vector<CigarElement> & dest, int pos, int start, int end, CigarElement elt);

    static bool isIndelAlignedTooFarLeft(Cigar* cigar, int leftmostAllowedAlignment);

    static bool cigarHasZeroSizeElement(Cigar* c);
};

#endif //MUTECT2CPP_MASTER_ALIGNMENTUTILS_H
