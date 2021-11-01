//
// Represents a pileup of reads at a given position
// Created by lhh on 10/29/21.
//

#ifndef MUTECT2CPP_MASTER_READPILEUP_H
#define MUTECT2CPP_MASTER_READPILEUP_H

#include "htslib/sam.h"
#include <vector>

class ReadPileup {
private:
    int tid;
    hts_pos_t pos;
    std::vector<bam1_t *> & pileupElements;

public:
    ReadPileup(int tid, hts_pos_t pos, std::vector<bam1_t *> & reads);

    hts_pos_t getPosition();

    std::vector<bam1_t* > & getPileupElements();
};


#endif //MUTECT2CPP_MASTER_READPILEUP_H
