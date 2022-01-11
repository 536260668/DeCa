//
// Created by 梦想家xixi on 2022/1/11.
//

#include "AlignmentContext.h"

AlignmentContext::AlignmentContext(std::vector<SAMRecord> &tumor, std::vector<SAMRecord> &normal, SimpleInterval &loc, int tid, SAMFileHeader* header) : tumor(std::move(tumor)), normal(std::move(normal)), loc(loc),
tid(tid), header(header){
}

int AlignmentContext::getReadNum() const{
    return normal.size() + tumor.size();
}

std::string AlignmentContext::getRefName() {
    return header->getSequenceDictionary().getSequences()[tid].getSequenceName();
}

int AlignmentContext::getPosition() const {
    return loc.getStart();
}

ReadPileup AlignmentContext::makeTumorPileup() {
    return {tid, loc.getStart(), tumor};
}

ReadPileup AlignmentContext::makeNormalPileup() {
    return {tid, loc.getStart(), normal};
}

bool AlignmentContext::isEmpty() const {
    return tumor.size() + normal.size() == 0;
}
