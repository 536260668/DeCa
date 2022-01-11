//
// Created by 梦想家xixi on 2022/1/11.
//

#ifndef MUTECT2CPP_MASTER_READPILEUP_H
#define MUTECT2CPP_MASTER_READPILEUP_H

#include "samtools/SAMRecord.h"

class ReadPileup {
private:
    int tid;
    int pos;
    std::vector<SAMRecord> reads;

public:
    ReadPileup(int tid, int pos, std::vector<SAMRecord>& reads);
    std::vector<SAMRecord> getPileupElements();
    int size() const {return reads.size();}
    int getPosition();
};


#endif //MUTECT2CPP_MASTER_READPILEUP_H
