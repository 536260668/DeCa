//
// Created by 梦想家xixi on 2022/1/11.
//

#include "ReadPileup.h"

ReadPileup::ReadPileup(int tid, int pos, std::vector<std::shared_ptr<SAMRecord>> &reads) : tid(tid), pos(pos), reads(reads){
}

std::vector<std::shared_ptr<SAMRecord>> ReadPileup::getPileupElements() {
    return reads;
}

int ReadPileup::getPosition() {
    return pos;
}
