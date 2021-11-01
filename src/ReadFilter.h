//
// The class used to filter out reads for Mutect2Cpp-master
// Created by lhh on 10/19/21.
//

#ifndef MUTECT2CPP_MASTER_READFILTER_H
#define MUTECT2CPP_MASTER_READFILTER_H

#include "htslib/sam.h"

class ReadFilter {  // TODO: add the filter left
public:
    static bool ReadLengthTest(bam1_t * read);
};


#endif //MUTECT2CPP_MASTER_READFILTER_H
