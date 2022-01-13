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
    std::vector<std::shared_ptr<SAMRecord>> tumor;
    std::vector<std::shared_ptr<SAMRecord>> normal;
    SimpleInterval loc;
    int tid;
    SAMFileHeader * header;

public:
    AlignmentContext(std::vector<std::shared_ptr<SAMRecord>> &tumor, std::vector<std::shared_ptr<SAMRecord>> &normal, SimpleInterval &loc, int tid, SAMFileHeader* header);
    AlignmentContext() {}
    int getReadNum() const;
    std::string getRefName();
    int getPosition() const;
    ReadPileup makeTumorPileup();
    ReadPileup makeNormalPileup();
    bool isEmpty() const;
    SimpleInterval& getLocation();
};


#endif //MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
