//
// Created by 梦想家xixi on 2021/12/18.
//

#ifndef MUTECT2CPP_MASTER_READUTILS_H
#define MUTECT2CPP_MASTER_READUTILS_H

#include "samtools/SAMRecord.h"

class ReadUtils {
public:
    static SAMRecord* emptyRead(SAMRecord* read);
    static void assertAttributeNameIsLegal(std::string& attributeName);
};


#endif //MUTECT2CPP_MASTER_READUTILS_H
