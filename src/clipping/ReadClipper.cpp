//
// Created by 梦想家xixi on 2021/12/18.
//

#include "ReadClipper.h"
#include "read/ReadUtils.h"

std::shared_ptr<SAMRecord> ReadClipper::hardClipToRegion(std::shared_ptr<SAMRecord> read, int refStart, int refStop) {
    int start = read->getStart();
    int stop = read->getEnd();
    return hardClipToRegion(read, refStart, refStop, start, stop);
}

std::shared_ptr<SAMRecord>
ReadClipper::hardClipToRegion(std::shared_ptr<SAMRecord> read, int refStart, int refStop, int alignmentStart, int alignmentStop) {
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

ReadClipper::ReadClipper(std::shared_ptr<SAMRecord> & read) : read(read), wasClipped(false){
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipBothEndsByReferenceCoordinates(const int left, const int right) {
    if(read->getLength() == 0 || left == right) {
        return ReadUtils::emptyRead(read);
    }
    std::shared_ptr<SAMRecord> leftTailRead = clipByReferenceCoordinates(right, -1, HARDCLIP_BASES, true);
    if(left > leftTailRead->getEnd()) {
        return ReadUtils::emptyRead(read);
    }
    ReadClipper clipper = ReadClipper(leftTailRead);
    return clipper.hardClipByReferenceCoordinatesLeftTail(left);
}

std::shared_ptr<SAMRecord>
ReadClipper::clipByReferenceCoordinates(int refStart, int refStop, ClippingRepresentation clippingOp, bool runAsserts) {
    if(read->getLength() == 0) {
        return read;
    }
    if(clippingOp == SOFTCLIP_BASES && read->isUnmapped()) {
        throw std::invalid_argument("Cannot soft-clip read by reference coordinates because it is unmapped");
    }
    int start;
    int stop;

    if(refStart < 0) {
        if(refStop < 0) {
            throw std::invalid_argument("Only one of refStart or refStop must be < 0");
        }
        start = 0;
        stop = ReadUtils::getReadCoordinateForReferenceCoordinate(read, refStop, LEFT_TAIL);
    } else {
        if(refStop >= 0) {
            throw std::invalid_argument("Either refStart or refStop must be < 0");
        }
        start = ReadUtils::getReadCoordinateForReferenceCoordinate(read, refStart, RIGHT_TAIL);
        stop = read->getLength() - 1;
    }

    if(start < 0 || stop > read->getLength() - 1) {
        throw std::invalid_argument("Trying to clip before the start or after the end of a read");
    }

    if(start > stop) {
        throw std::invalid_argument("START > STOP -- this should never happen, please check read");
    }

    if(start > 0 && stop < read->getLength() - 1) {
        throw std::invalid_argument("Trying to clip the middle of the read");
    }
    addOp(ClippingOp(start, stop));
    std::shared_ptr<SAMRecord> clippedRead = clipRead(clippingOp, runAsserts);
    ops.clear();
    return clippedRead;
}

void ReadClipper::addOp(const ClippingOp &op) {
    ops.emplace_back(op);
}

std::shared_ptr<SAMRecord> ReadClipper::clipRead(ClippingRepresentation algorithm, bool runAsserts) {
    Mutect2Utils::validateArg(algorithm != NULL_ClippingRepresentation, "null is not allowed there.");
    if(ops.empty()) {
        return read;
    }
    std::shared_ptr<SAMRecord> clippedRead = read;
    for(ClippingOp op : ops) {
        int readLength = clippedRead->getLength();
        if(op.start < readLength) {
            ClippingOp fixedOperation = op;
            if(op.stop >= readLength) {
                fixedOperation = ClippingOp(op.start, readLength - 1);
            }
            clippedRead = fixedOperation.apply(algorithm, clippedRead, runAsserts);
        }
    }
    wasClipped = true;
    ops.clear();
    if(clippedRead->isEmpty()) {
        return ReadUtils::emptyRead(clippedRead);
    }
    return clippedRead;
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipByReferenceCoordinatesLeftTail(int refStop) {
    return clipByReferenceCoordinates(-1, refStop, HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipBothEndsByReferenceCoordinates(std::shared_ptr<SAMRecord> read, int left, int right) {
    return ReadClipper(read).hardClipBothEndsByReferenceCoordinates(left, right);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipByReferenceCoordinatesLeftTail(std::shared_ptr<SAMRecord> read, int refStop) {
    return ReadClipper(read).clipByReferenceCoordinates(-1, refStop, HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipByReferenceCoordinatesRightTail(std::shared_ptr<SAMRecord> read, int refStart) {
    return ReadClipper(read).clipByReferenceCoordinates(refStart, -1, HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::clipLowQualEnds(ClippingRepresentation algorithm, uint8_t lowQual) {
    if(read->isEmpty())
        return read;

    int readLength = read->getLength();
    int leftClipIndex = 0;
    int rightClipIndex = readLength - 1;

    while (rightClipIndex >= 0 && read->getBaseQuality(rightClipIndex) <= lowQual) {
        rightClipIndex--;
    }
    while (leftClipIndex < readLength && read->getBaseQuality(leftClipIndex) <= lowQual) {
        leftClipIndex++;
    }
    if (leftClipIndex > rightClipIndex) {
        return ReadUtils::emptyRead(read);
    }

    if(rightClipIndex < readLength - 1) {
        addOp(ClippingOp(rightClipIndex + 1, readLength - 1));
    }

    if (leftClipIndex > 0 ) {
        addOp(ClippingOp(0, leftClipIndex - 1));
    }
    return clipRead(algorithm, true);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipLowQualEnds(uint8_t lowQual) {
    return clipLowQualEnds(HARDCLIP_BASES, lowQual);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipLowQualEnds(std::shared_ptr<SAMRecord> read, uint8_t lowQual) {
    return ReadClipper(read).hardClipLowQualEnds(lowQual);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipSoftClippedBases(std::shared_ptr<SAMRecord> read) {
    return ReadClipper(read).hardClipSoftClippedBases();
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipSoftClippedBases() {
    if(read->isEmpty()) {
        return read;
    }

    int readIndex = 0;
    int cutLeft = -1;            // first position to hard clip (inclusive)
    int cutRight = -1;           // first position to hard clip (inclusive)
    bool rightTail = false;

    for(CigarElement cigarElement : read->getCigarElements()) {
        if (cigarElement.getOperator() == S) {
            if (rightTail) {
                cutRight = readIndex;
            }
            else {
                cutLeft = readIndex + cigarElement.getLength() - 1;
            }
        }
        else if (cigarElement.getOperator() != H) {
            rightTail = true;
        }

        if (CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator()) ) {
            readIndex += cigarElement.getLength();
        }
    }

    if (cutRight >= 0) {
        addOp(ClippingOp(cutRight, read->getLength() - 1));
    }
    if (cutLeft >= 0) {
        addOp(ClippingOp(0, cutLeft));
    }
    return clipRead(HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::revertSoftClippedBases(std::shared_ptr<SAMRecord> read) {
    return ReadClipper(read).revertSoftClippedBases();
}

std::shared_ptr<SAMRecord> ReadClipper::revertSoftClippedBases() {
    if (read->isEmpty()) {
        return read;
    }
    addOp(ClippingOp(0, 0));
    return clipRead(REVERT_SOFTCLIPPED_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipAdaptorSequence(std::shared_ptr<SAMRecord> read) {
    return ReadClipper(read).hardClipAdaptorSequence();
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipAdaptorSequence() {
    int adaptorBoundary = read->getAdaptorBoundary();

    if (adaptorBoundary == ReadUtils::CANNOT_COMPUTE_ADAPTOR_BOUNDARY || !ReadUtils::isInsideRead(read, adaptorBoundary)) {
        return read;
    }

    return read->isReverseStrand() ? hardClipByReferenceCoordinatesLeftTail(adaptorBoundary) : hardClipByReferenceCoordinatesRightTail(read, adaptorBoundary);
}
