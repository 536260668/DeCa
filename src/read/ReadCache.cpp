//
// Created by 梦想家xixi on 2022/1/11.
//

#include "ReadCache.h"
#include "iostream"
#include "ReadUtils.h"

ReadCache::ReadCache(aux_t **data, std::vector<char*> & bam_name) : data(data), tid(0), bam_name(bam_name){
    bam1_t * b;
    b = bam_init1();
    std::string region = data[0]->header->getSequenceDictionary().getSequences()[0].getSequenceName();
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        int count = 0;
        idx = sam_index_load(data[i]->fp, bam_name[i]);
        if(idx == 0)
            throw std::invalid_argument("random alignment retrieval only works for indexed BAM or CRAM files.");
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            SAMRecord read(b, data[i]->header);
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0)
                    normalReads.emplace_back(read);
                else
                    tumorReads.emplace_back(read);
                count++;
            }
            if(count > 500)
                break;
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);
    start = end = currentPose = 0;
}

ReadCache::ReadCache(aux_t **data, std::vector<char *> &bam_name, int tid, const std::string& region) : tid(tid), data(data), bam_name(bam_name){
    bam1_t * b;
    b = bam_init1();
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        idx = sam_index_load(data[i]->fp, bam_name[i]);
        if(idx == 0)
            throw std::invalid_argument("random alignment retrieval only works for indexed BAM or CRAM files.");
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            SAMRecord read(b, data[i]->header);
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0)
                    normalReads.emplace_back(read);
                else
                    tumorReads.emplace_back(read);
            }
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);
    int i = region.find(':') + 1;
    start = end = 0;
    while(region[i] >= '0' && region[i] <= '9') {
        start = start * 10 + region[i] - '0';
        i++;
    }
    currentPose = start - 1;
    i = region.find('-') + 1;
    while(region[i] >= '0' && region[i] <= '9') {
        end = end * 10 + region[i] - '0';
        i++;
    }
}

ReadCache::~ReadCache() {
    hts_idx_destroy(idx);
}

int ReadCache::getNextPos() {
    int nextPose = currentPose + 1;
    if(nextPose > data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceLength()) {
        throw std::invalid_argument("please check first");
    }
    while(nextPose > end || (tumorReads.empty() && normalReads.empty())) {
        advanceLoad();
    }
    if((tumorReads.empty() || tumorReads.front().getStart() > nextPose) && (normalReads.empty() || normalReads.front().getStart() > nextPose)){
        if(tumorReads.empty()){
            nextPose = normalReads.front().getStart();
        } else if (normalReads.empty()) {
            nextPose = tumorReads.front().getStart();
        } else
            nextPose = std::min(tumorReads.front().getStart(), normalReads.front().getStart());
    }
    currentPose = nextPose;
    return nextPose;
}

bool ReadCache::hasNextPos() {
    int nextPose = currentPose + 1;
    return nextPose <= data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceLength();
}

void ReadCache::advanceLoad() {
    tumorReads.clear();
    normalReads.clear();
    start = end + 1;
    end = start + 100000 -1;
    std::string region = data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' +
            std::to_string(start) + '-' + std::to_string(end);
    bam1_t * b;
    b = bam_init1();
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            SAMRecord read(b, data[i]->header);
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0)
                    normalReads.emplace_back(read);
                else
                    tumorReads.emplace_back(read);
            }
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);
}

AlignmentContext ReadCache::getAlignmentContext() {
    getNextPos();
    std::vector<SAMRecord> tumor;
    std::vector<SAMRecord> normal;
    while(tumorReads.front().getEnd() < currentPose)
        tumorReads.pop_front();
    while(normalReads.front().getEnd() < currentPose)
        normalReads.pop_front();
    std::deque<SAMRecord>::iterator iter = tumorReads.begin();
    while(iter != tumorReads.end() && iter->getStart() < currentPose) {
        if(iter->getEnd() >= currentPose && !ReadUtils::isBaseInsideAdaptor(&(*iter), currentPose)) {
            tumor.emplace_back(*iter);
        }
        iter++;
    }
    iter = normalReads.begin();
    while(iter != normalReads.end() && iter->getStart() < currentPose) {
        if(iter->getEnd() > currentPose && !ReadUtils::isBaseInsideAdaptor(&(*iter), currentPose)) {
            normal.emplace_back(*iter);
        }
        iter++;
    }
    SimpleInterval loc(data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceName(), currentPose, currentPose);
    return {tumor, normal, loc, tid, data[0]->header};
}
