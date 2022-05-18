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
#include "transformers/PalindromeArtifactClipReadTransformer.h"
#include "transformers/BQSRReadTransformer.h"
#include "pileRead.h"

#define REGION_SIZE 1000000

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    sam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    SAMFileHeader * header;
} aux_t;


class AlignmentContext;

class ReadCache {
private:
    std::queue<pileRead*> tumorReads;
    std::queue<pileRead*> normalReads;
    std::list<std::shared_ptr<SAMRecord>> tumorReadsForRegion;
    std::list<std::shared_ptr<SAMRecord>> normalReadsForRegion;
    std::list<pileRead*> tumorReadsForAlignment;
    std::list<pileRead*> normalReadsForAlignment;
    std::list<pileRead*> tumorCache;
    std::list<pileRead*> normalCache;
    std::vector<char*> bam_name;
    aux_t ** data;
    int tid;
    int ExtendedStart;
    int ExtendedEnd;
    int chr_len;    // the total length of current chromosome
    std::string chr_name;   // the name of current chromosome
    std::vector<hts_idx_t *> hts_idxes;
    int currentPose;
    PalindromeArtifactClipReadTransformer readTransformer;
    BQSRReadTransformer * tumorTransformer;
    BQSRReadTransformer * normalTransformer;
    bool bqsr_within_mutect;

    // read data in a specific region and add it to the cache
    void readData(const string & region);

    void advanceLoad();

    // clear all the reads in the cache
    void clear();

public:

    ReadCache(aux_t** data, std::vector<char*> & bam_name, std::shared_ptr<ReferenceCache> & cache);
    ReadCache(aux_t **data, std::vector<char *> &bam_name, int tid, int start, int end, int maxAssemblyRegionSize, std::shared_ptr<ReferenceCache> & cache, bool bqsr_within_mutect = false, BQSRReadTransformer * tumorTransformer = nullptr, BQSRReadTransformer * normalTransformer = nullptr);
    int getNextPos();
    bool hasNextPos();
    void InsertPileToAlignment(pileRead* stopPos, std::list<pileRead*> &);
    bool InsertPileToCache(pileRead* stopPos, std::list<pileRead*> & toAdd);
    AlignmentContext getAlignmentContext();
    std::vector<std::shared_ptr<SAMRecord>> getReadsForRegion(AssemblyRegion & region);
    static pileRead * getpileRead(const std::shared_ptr<SAMRecord> & read);
    ~ReadCache();
};


#endif //MUTECT2CPP_MASTER_READCACHE_H
