//
// Created by 梦想家xixi on 2022/1/11.
//

#include "ReadPileup.h"

ReadPileup::ReadPileup(int tid, int pos, std::vector<SAMRecord> &reads) : tid(tid), pos(pos), reads(reads){
}

std::vector<SAMRecord> ReadPileup::getPileupElements() {
    return reads;
}

int ReadPileup::getPosition() {
    return pos;
}
