//
// Created by 梦想家xixi on 2022/1/7.
//

#ifndef MUTECT2CPP_MASTER_LOCUSITERATORBYSTATE_H
#define MUTECT2CPP_MASTER_LOCUSITERATORBYSTATE_H

#include "AlignmentContext.h"
#include "samtools/SAMRecord.h"

class LocusIteratorByState {
private:
    int n;
    SAMFileHeader* header;
    std::vector<sam_hdr_t *> & headers;
    int *n_plp;
    int pos;
    int tid;
    bam_pileup1_t **plp;
    std::vector<SAMRecord> tumor;
    std::vector<SAMRecord> normal;

public:
    LocusIteratorByState(std::vector<sam_hdr_t *> & headers);
    LocusIteratorByState(int n, int * nplp, bam_pileup1_t** plp, SAMFileHeader* header, std::vector<sam_hdr_t *> & headers, int pos, int tid);
    ~LocusIteratorByState();
    AlignmentContext loadAlignmentContext();
};


#endif //MUTECT2CPP_MASTER_LOCUSITERATORBYSTATE_H
