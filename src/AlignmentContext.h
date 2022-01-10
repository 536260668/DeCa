//
// Bundles together a pileup and a location
// Different from GATK4, this class is implemented without using pileup
//
// Created by lhh on 10/22/21.
//

#ifndef MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
#define MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H

#include <set>
#include "htslib/sam.h"
#include "ReadPileup.h"

class AlignmentContext {
private:
    std::vector<sam_hdr_t *> & headers;
    int tid;
    hts_pos_t pos;
    int n;
    int * n_plp;
    bam_pileup1_t ** plp;
    bool flag;

public:
    AlignmentContext(std::vector<sam_hdr_t *> & headers, int tid, hts_pos_t pos, int n, int * n_plp, bam_pileup1_t ** plp, bool flag = false);

    ~AlignmentContext();

    hts_pos_t getPosition() const;

    int getTid() const;

    int getReadNum() const;

    const char * getRefName();

    // added by lhh, used to get the tumor pileup from the reads
    ReadPileup makeTumorPileup(std::set<std::string> & normalSamples);

    ReadPileup makeNormalPileup(std::set<std::string> & normalSamples);

    bool isEmpty();
};


#endif //MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
