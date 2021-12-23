//
// Created by 梦想家xixi on 2021/12/18.
//

#ifndef MUTECT2CPP_MASTER_READUTILS_H
#define MUTECT2CPP_MASTER_READUTILS_H

#include "samtools/SAMRecord.h"
#include "samtools/SAMFileHeader.h"

enum ClippingTail {
    LEFT_TAIL,
    RIGHT_TAIL
};

class ReadUtils {
public:
    static SAMRecord* emptyRead(SAMRecord* read);
    static void assertAttributeNameIsLegal(std::string& attributeName);
    static int getReadCoordinateForReferenceCoordinate(SAMRecord* read, int refCoord, ClippingTail tail);
    static int getReadCoordinateForReferenceCoordinate(int alignmentStart, Cigar* cigar, int refCoord, ClippingTail tail, bool allowGoalNotReached);
    static CigarElement* readStartsWithInsertion(Cigar* cigarForRead);
    static CigarElement* readStartsWithInsertion(Cigar* cigarForRead, bool ignoreSoftClipOps);
    static int getSoftStart(SAMRecord* read);
    static int getSoftEnd(SAMRecord* read);
    static bool hasBaseIndelQualities(SAMRecord* read);
    static std::string BQSR_BASE_INSERTION_QUALITIES;
    static std::string BQSR_BASE_DELETION_QUALITIES;
    static uint8_t * getExistingBaseInsertionQualities(SAMRecord* read, int & length);
    static uint8_t * getExistingBaseDeletionQualities(SAMRecord* read, int & length);
    static uint8_t * getBaseInsertionQualities(SAMRecord* read, int & length);
    static uint8_t * getBaseDeletionQualities(SAMRecord* read, int & length);
    static void setInsertionBaseQualities(SAMRecord* read, uint8_t* quals, int length);
    static void setDeletionBaseQualities(SAMRecord* read, uint8_t* quals, int length);
    static int getAssignedReferenceIndex(SAMRecord* read, SAMFileHeader* header);
    static int getSAMFlagsForRead(SAMRecord* read);
    static int getMateReferenceIndex(SAMRecord* read, SAMFileHeader* header);

private:
    static std::pair<int, bool> getReadCoordinateForReferenceCoordinate(int alignmentStart, Cigar* cigar, int refCoord, bool allowGoalNotReached);
    static const int CLIPPING_GOAL_NOT_REACHED = -1;
};


#endif //MUTECT2CPP_MASTER_READUTILS_H
