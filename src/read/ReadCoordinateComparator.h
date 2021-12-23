//
// Created by 梦想家xixi on 2021/12/23.
//

#ifndef MUTECT2CPP_MASTER_READCOORDINATECOMPARATOR_H
#define MUTECT2CPP_MASTER_READCOORDINATECOMPARATOR_H


#include "samtools/SAMFileHeader.h"
#include "samtools/SAMRecord.h"

class ReadCoordinateComparator {
private:
    SAMFileHeader* header;

public:
    explicit ReadCoordinateComparator(SAMFileHeader* header);
    int compare(SAMRecord* first, SAMRecord* second);
    static int compareCoordinates(SAMRecord* first, SAMRecord* second, SAMFileHeader* header);
};


#endif //MUTECT2CPP_MASTER_READCOORDINATECOMPARATOR_H
