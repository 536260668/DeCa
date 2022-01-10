//
// Created by lhh on 10/29/21.
//

#include "ReadPileup.h"

ReadPileup::ReadPileup(int tid, hts_pos_t pos, std::vector<bam1_t *> &reads): tid(tid), pos(pos), pileupElements(reads)
{}

hts_pos_t ReadPileup::getPosition()
{
    return pos;
}

std::vector<bam1_t* > & ReadPileup::getPileupElements()
{
    return pileupElements;
}

int ReadPileup::size() {
    return pileupElements.size();
}
