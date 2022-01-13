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
    std::shared_ptr<SAMRecord> read;
    bool wasClipped;
    std::vector<ClippingOp> ops;
    explicit ReadClipper(std::shared_ptr<SAMRecord> & read);
    void addOp(const ClippingOp & op);
    static std::shared_ptr<SAMRecord>hardClipToRegion(std::shared_ptr<SAMRecord> read, int refStart, int refStop);
    static std::shared_ptr<SAMRecord> hardClipBothEndsByReferenceCoordinates(std::shared_ptr<SAMRecord> read, int left, int right);
    static std::shared_ptr<SAMRecord> hardClipByReferenceCoordinatesLeftTail(std::shared_ptr<SAMRecord> read, int refStop);
    static std::shared_ptr<SAMRecord> hardClipByReferenceCoordinatesRightTail(std::shared_ptr<SAMRecord> read, int refStop);

private:
    static std::shared_ptr<SAMRecord> hardClipToRegion(std::shared_ptr<SAMRecord> read, int refStart, int refStop, int alignmentStart, int alignmentStop);
    std::shared_ptr<SAMRecord> hardClipBothEndsByReferenceCoordinates(int left, int right);
    std::shared_ptr<SAMRecord> clipByReferenceCoordinates(int refStart, int refStop, ClippingRepresentation clippingOp, bool runAsserts);
    std::shared_ptr<SAMRecord> clipRead(ClippingRepresentation algorithm, bool runAsserts);
    std::shared_ptr<SAMRecord> hardClipByReferenceCoordinatesLeftTail(int refStop);
};


#endif //MUTECT2CPP_MASTER_READCLIPPER_H
