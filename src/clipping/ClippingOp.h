//
// Created by 梦想家xixi on 2021/12/21.
//

#ifndef MUTECT2CPP_MASTER_CLIPPINGOP_H
#define MUTECT2CPP_MASTER_CLIPPINGOP_H

#include "samtools/SAMRecord.h"
#include "ReadClipper.h"
#include "ClippingRepresentation.h"


class ClippingOp {
public:
    int start;
    int stop;
    ClippingOp(int start, int stop);
    SAMRecord* apply(ClippingRepresentation algorithm, SAMRecord* originalRead, bool runAsserts);

private:
    SAMRecord* applyHardClipBases(SAMRecord* read, int start, int stop);
    class CigarShift {
    public:
        CigarShift(Cigar* cigar, int shiftFromStart, int shiftFromEnd);

        int shiftFromStart;
        int shiftFromEnd;
        Cigar* cigar;
    };
    CigarShift* hardClipCigar(Cigar* cigar, int start, int stop);
    int calculateHardClippingAlignmentShift(CigarElement & cigarElement, int clippedLength);
    CigarShift* cleanHardClippedCigar(Cigar* cigar);
    static int calculateAlignmentStartShift(Cigar* oldCigar, Cigar* newCigar);
    static int calcHardSoftOffset(Cigar* cigar);
};



#endif //MUTECT2CPP_MASTER_CLIPPINGOP_H
