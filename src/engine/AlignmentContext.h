//
// Created by 梦想家xixi on 2022/1/11.
//

#ifndef MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
#define MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H

#include "samtools/SAMRecord.h"
#include "utils/ReadPileup.h"
#include "SimpleInterval.h"

class AlignmentContext {
private:
    std::vector<SAMRecord> tumor;
    std::vector<SAMRecord> normal;
    SimpleInterval loc;
    int tid;
    SAMFileHeader * header;

public:
    AlignmentContext(std::vector<SAMRecord> &tumor, std::vector<SAMRecord> &normal, SimpleInterval &loc, int tid, SAMFileHeader* header);
    int getReadNum() const;
    std::string getRefName();
    int getPosition() const;
    ReadPileup makeTumorPileup();
    ReadPileup makeNormalPileup();
    bool isEmpty() const;
};


#endif //MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
