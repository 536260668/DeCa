//
// Created by 梦想家xixi on 2021/12/18.
//

#ifndef MUTECT2CPP_MASTER_READCLIPPER_H
#define MUTECT2CPP_MASTER_READCLIPPER_H

#include "samtools/SAMRecord.h"
#include "ClippingOp.h"
#include "ClippingRepresentation.h"

class ClippingOp;

class ReadClipper {
public:
    SAMRecord* read;
    bool wasClipped;
    std::vector<ClippingOp> ops;
    ReadClipper(SAMRecord* read);
    void addOp(const ClippingOp & op);
    static SAMRecord* hardClipToRegion(SAMRecord * read, int refStart, int refStop);
    static SAMRecord* hardClipBothEndsByReferenceCoordinates(SAMRecord* read, int left, int right);
    static SAMRecord* hardClipByReferenceCoordinatesLeftTail(SAMRecord* read, int refStop);
    static SAMRecord* hardClipByReferenceCoordinatesRightTail(SAMRecord* read, int refStop);

private:
    static SAMRecord* hardClipToRegion(SAMRecord * read, int refStart, int refStop, int alignmentStart, int alignmentStop);
    SAMRecord* hardClipBothEndsByReferenceCoordinates(int left, int right);
    SAMRecord* clipByReferenceCoordinates(int refStart, int refStop, ClippingRepresentation clippingOp, bool runAsserts);
    SAMRecord* clipRead(ClippingRepresentation algorithm, bool runAsserts);
    SAMRecord* hardClipByReferenceCoordinatesLeftTail(int refStop);
};


#endif //MUTECT2CPP_MASTER_READCLIPPER_H
