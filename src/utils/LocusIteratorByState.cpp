//
// Created by 梦想家xixi on 2022/1/7.
//

#include "LocusIteratorByState.h"
#include "read/ReadUtils.h"
#include <iostream>

LocusIteratorByState::LocusIteratorByState(std::vector<sam_hdr_t *> & headers) : n(0), n_plp(nullptr), plp(nullptr), header(nullptr), pos(0), headers(headers), tid(-1){

}

LocusIteratorByState::LocusIteratorByState(int n, int *nplp, bam_pileup1_t **plp, SAMFileHeader *header, std::vector<sam_hdr_t *> & headers, int pos, int tid) : n(n), header(header), n_plp(nplp), plp(plp), pos(pos), headers(headers),
tid(tid){
    Mutect2Utils::validateArg(n == 2, "input error");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < nplp[i]; j++) {
            //std::cout << nplp[i] << std::endl;
            if(i == 0) {
                normal.emplace_back(SAMRecord(plp[i][j].b, header, false));
            } else {
                tumor.emplace_back(SAMRecord(plp[i][j].b, header, false));
            }
        }
    }
}

AlignmentContext LocusIteratorByState::loadAlignmentContext() {
    int * new_n_plp = new int[n]{0};
    bam_pileup1_t ** new_plp = new bam_pileup1_t*[n];
    std::vector<SAMRecord> new_normal;
    std::vector<SAMRecord> new_tumor;
    std::vector<bam_pileup1_t> tmp;
    for(int j = 0; j < n_plp[0]; j++) {
        if(!ReadUtils::isBaseInsideAdaptor(&normal[j], pos)) {
            tmp.emplace_back(plp[0][j]);
            new_n_plp[0]++;
        }
    }
    bam_pileup1_t * afterFilter = new bam_pileup1_t[new_n_plp[0]];
    for(int i = 0; i < new_n_plp[0]; i++){
        afterFilter[i] = tmp[i];
    }
    new_plp[0] = afterFilter;
    tmp.clear();
    for(int j = 0; j < n_plp[1]; j++) {
        if(!ReadUtils::isBaseInsideAdaptor(&tumor[j], pos)) {
            tmp.emplace_back(plp[1][j]);
            new_n_plp[1]++;
        }
    }
    afterFilter = new bam_pileup1_t[new_n_plp[1]];
    for(int i = 0; i < new_n_plp[1]; i++){
        afterFilter[i] = tmp[i];
    }
    new_plp[1] = afterFilter;
//    AlignmentContext ret(headers, tid, pos, n, new_n_plp, new_plp, true);
    return {headers, tid, pos, n, new_n_plp, new_plp, true};
}

LocusIteratorByState::~LocusIteratorByState() {
//    tumor.clear();
}




