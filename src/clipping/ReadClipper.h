//
// Created by 梦想家xixi on 2021/12/18.
//

#ifndef MUTECT2CPP_MASTER_READCLIPPER_H
#define MUTECT2CPP_MASTER_READCLIPPER_H

#include "samtools/SAMRecord.h"

class ReadClipper {
public:
    static SAMRecord* hardClipToRegion(SAMRecord * read, int refStart, int refStop);
    SAMRecord* read;
    bool wasClipped;

private:
    static SAMRecord* hardClipToRegion(SAMRecord * read, int refStart, int refStop, int alignmentStart, int alignmentStop);
};


#endif //MUTECT2CPP_MASTER_READCLIPPER_H
