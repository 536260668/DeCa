//
// Created by 梦想家xixi on 2021/12/16.
//

#ifndef MUTECT2CPP_MASTER_SAMSEQUENCEDICTIONARY_H
#define MUTECT2CPP_MASTER_SAMSEQUENCEDICTIONARY_H

#include "SAMSequenceRecord.h"
#include <map>
#include <vector>

class SAMSequenceDictionary {
private:
    std::vector<SAMSequenceRecord> mSequences;
    std::map<std::string, SAMSequenceRecord> mSequenceMap;

public:
    SAMSequenceDictionary() = default;
    void addSequence(SAMSequenceRecord sequenceRecord);
    SAMSequenceRecord & getSequence(const std::string& name);
    int getSequenceIndex(std::string & sequenceName);
};


#endif //MUTECT2CPP_MASTER_SAMSEQUENCEDICTIONARY_H
