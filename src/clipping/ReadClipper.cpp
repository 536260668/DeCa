//
// Created by 梦想家xixi on 2021/12/18.
//

#include "ReadClipper.h"
#include "read/ReadUtils.h"

SAMRecord *ReadClipper::hardClipToRegion(SAMRecord *read, int refStart, int refStop) {
    int start = read->getStart();
    int stop = read->getEnd();
    return hardClipToRegion(read, refStart, refStop, start, stop);
}

SAMRecord *
ReadClipper::hardClipToRegion(SAMRecord *read, int refStart, int refStop, int alignmentStart, int alignmentStop) {
    if(alignmentStart <= refStop && alignmentStop >= refStart) {
        if(alignmentStart < refStart && alignmentStop > refStop) {
           return hardClipBothEndsByReferenceCoordinates(read, refStart - 1, refStop + 1);
        } else if (alignmentStart < refStart) {
            return hardClipByReferenceCoordinatesLeftTail(read, refStart - 1);
        } else if (alignmentStop > refStop) {
            return hardClipByReferenceCoordinatesRightTail(read, refStop + 1);
        }
        return read;
    } else {
        return ReadUtils::emptyRead(read);
    }
}
