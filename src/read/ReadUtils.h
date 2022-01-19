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
    static const int CANNOT_COMPUTE_ADAPTOR_BOUNDARY = INT32_MIN;
    static std::shared_ptr<SAMRecord> emptyRead(std::shared_ptr<SAMRecord> & read);
    static void assertAttributeNameIsLegal(std::string& attributeName);
    static int getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> & read, int refCoord, ClippingTail tail);
    static int getReadCoordinateForReferenceCoordinate(int alignmentStart, Cigar* cigar, int refCoord, ClippingTail tail, bool allowGoalNotReached);
    static CigarElement* readStartsWithInsertion(Cigar* cigarForRead);
    static CigarElement* readStartsWithInsertion(Cigar* cigarForRead, bool ignoreSoftClipOps);
    static int getSoftStart(std::shared_ptr<SAMRecord> & read);
    static int getSoftEnd(std::shared_ptr<SAMRecord> & read);
    static bool hasBaseIndelQualities(std::shared_ptr<SAMRecord> & read);
    static std::string BQSR_BASE_INSERTION_QUALITIES;
    static std::string BQSR_BASE_DELETION_QUALITIES;
    static uint8_t * getExistingBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static uint8_t * getExistingBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static uint8_t * getBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static uint8_t * getBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static void setInsertionBaseQualities(std::shared_ptr<SAMRecord> & read, uint8_t* quals, int length);
    static void setDeletionBaseQualities(std::shared_ptr<SAMRecord> & read, uint8_t* quals, int length);
    static int getAssignedReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader* header);
    static int getSAMFlagsForRead(std::shared_ptr<SAMRecord> & read);
    static int getMateReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader* header);
    static bool alignmentAgreesWithHeader(SAMFileHeader* header, std::shared_ptr<SAMRecord> & read);
    static int getReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader* header);
    static bool hasWellDefinedFragmentSize(std::shared_ptr<SAMRecord> & read);
    static int getAdaptorBoundary(std::shared_ptr<SAMRecord> & read);
    static bool isBaseInsideAdaptor(std::shared_ptr<SAMRecord> & read, long basePos);
    static bool isInsideRead(std::shared_ptr<SAMRecord> & read, int referenceCoordinate);

private:
    static std::pair<int, bool> getReadCoordinateForReferenceCoordinate(int alignmentStart, Cigar* cigar, int refCoord, bool allowGoalNotReached);
    static const int CLIPPING_GOAL_NOT_REACHED;
};


#endif //MUTECT2CPP_MASTER_READUTILS_H
