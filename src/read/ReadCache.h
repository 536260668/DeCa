//
// Created by 梦想家xixi on 2022/1/11.
//

#ifndef MUTECT2CPP_MASTER_READCACHE_H
#define MUTECT2CPP_MASTER_READCACHE_H

#include "htslib/sam.h"
#include "samtools/SAMFileHeader.h"
#include "samtools/SAMRecord.h"
#include <deque>
#include "engine/AlignmentContext.h"
#include "ReadFilter.h"

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    sam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ; // mapQ filter;
    uint32_t flags;// read filtering flags
    SAMFileHeader * header;
} aux_t;

class ReadCache {
private:
    std::deque<SAMRecord> tumorReads;
    std::deque<SAMRecord> normalReads;
    std::vector<char*> bam_name;
    aux_t ** data;
    int tid;
    int start;
    int end;
    hts_idx_t *idx;
    int currentPose;
    void advanceLoad();

public:
    ReadCache(aux_t** data, std::vector<char*> & bam_name);
    ReadCache(aux_t** data, std::vector<char*> & bam_name, int tid, const std::string&);
    int getNextPos();
    bool hasNextPos();
    AlignmentContext getAlignmentContext();

    ~ReadCache();
};


#endif //MUTECT2CPP_MASTER_READCACHE_H
