//
// Created by lhh on 10/22/21.
//

#include "assert.h"
#include "AlignmentContext.h"
#include <iostream>

AlignmentContext::AlignmentContext(std::vector<sam_hdr_t *> & headers, int tid, hts_pos_t pos, int n, int *n_plp, bam_pileup1_t **plp, bool flag) : headers(headers),
                                                                                                                                                    tid(tid), pos(pos), n(n), n_plp(n_plp), plp(plp), flag(flag)
{

}

hts_pos_t AlignmentContext::getPosition() const
{
    return pos;
}

int AlignmentContext::getTid() const
{
    return tid;
}

int AlignmentContext::getReadNum() const
{
    int num = 0;
    for(int i=0; i<n; i++)
    {
        num += n_plp[i];
    }
    return num;
}

const char *AlignmentContext::getRefName()
{
    return sam_hdr_tid2name(headers[0], tid);   //---any problem?
}

ReadPileup AlignmentContext::makeTumorPileup(std::set<std::string> & normalSamples)
{
    //---Maybe this method can be refactored
    std::vector<bam1_t *> reads;
    bool IsTumor = true;
    kstring_t pu_val = KS_INITIALIZE;
    for(int i=0; i<n; i++)
    {
        IsTumor = true;
        int res = sam_hdr_find_tag_id(headers[i], "RG", NULL, NULL, "SM", &pu_val);
        assert(!res);
        for(const std::string & sample : normalSamples)
        {
            if(!strcmp(sample.c_str(), pu_val.s))
            {
                IsTumor = false;
            }
        }
        if(IsTumor)
        {
            for(int j=0; j<n_plp[i]; j++)
            {
                reads.emplace_back(plp[i][j].b);
            }
        }
    }
    return ReadPileup(tid, pos, reads);
}

bool AlignmentContext::isEmpty() {
    int count = 0;
    for(int i = 0; i < n; i++) {
        count += n_plp[i];
    }
    return count == 0;
}

AlignmentContext::~AlignmentContext() {
    if(flag) {
        for(int i = 0; i < n; i++) {
            if(n_plp[i] > 0)
                delete[] plp[i];
        }
        delete[] plp;
        delete[] n_plp;
    }
}

ReadPileup AlignmentContext::makeNormalPileup(std::set<std::string> &normalSamples) {
    std::vector<bam1_t *> reads;
    bool IsNormal = true;
    kstring_t pu_val = KS_INITIALIZE;
    for(int i=0; i<n; i++)
    {
        IsNormal = true;
        int res = sam_hdr_find_tag_id(headers[i], "RG", NULL, NULL, "SM", &pu_val);
        assert(!res);
        for(const std::string & sample : normalSamples)
        {
            if(strcmp(sample.c_str(), pu_val.s))
            {
                IsNormal = false;
            }
        }
        if(IsNormal)
        {
            for(int j=0; j<n_plp[i]; j++)
            {
                reads.emplace_back(plp[i]->b);
            }
        }
    }
    return {tid, pos, reads};
}
