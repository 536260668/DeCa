//
// Created by 梦想家xixi on 2022/1/11.
//

#ifndef MUTECT2CPP_MASTER_READCACHE_H
#define MUTECT2CPP_MASTER_READCACHE_H

#include "htslib/sam.h"
#include "samtools/SAMFileHeader.h"
#include "samtools/SAMRecord.h"
#include <list>
#include <queue>
#include "engine/AlignmentContext.h"
#include "ReadFilter.h"
#include "AssemblyRegion.h"

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
    std::queue<std::shared_ptr<SAMRecord>> tumorReads;
    std::queue<std::shared_ptr<SAMRecord>> normalReads;
    std::list<std::shared_ptr<SAMRecord>> tumorReadsForRegion;
    std::list<std::shared_ptr<SAMRecord>> normalReadsForRegion;
    std::list<std::shared_ptr<SAMRecord>> tumorReadsForAlignment;
    std::list<std::shared_ptr<SAMRecord>> normalReadsForAlignment;
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
    std::vector<std::shared_ptr<SAMRecord>> getReadsForRegion(AssemblyRegion & region);
    ~ReadCache();
};


#endif //MUTECT2CPP_MASTER_READCACHE_H
