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
            std::shared_ptr<SAMRecord> read(new SAMRecord(b, data[i]->header));
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(read);
                }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(read);
                }
                count++;
            }
            if(count > 500)
                break;
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);
    hts_idx_destroy(idx);
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
            std::shared_ptr<SAMRecord> read(new SAMRecord(b, data[i]->header));
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(read);
                    }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(read);
                }
            }
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);
    hts_idx_destroy(idx);
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

}

int ReadCache::getNextPos() {
    int nextPose = currentPose + 1;
    if(nextPose > data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceLength()) {
        throw std::invalid_argument("please check first");
    }
    while(nextPose > end) {
        advanceLoad();
    }
    currentPose = nextPose;
    return nextPose;
}

bool ReadCache::hasNextPos() {
    int nextPose = currentPose + 1;
    return nextPose <= data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceLength();
}

void ReadCache::advanceLoad() {
    if(!normalReads.empty() && !tumorReads.empty())
        throw std::logic_error("error");
    start = end + 1;
    end = start + 1000000 -1;

    std::string region = data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' +
            std::to_string(start+1) + '-' + std::to_string(end);
    bam1_t * b;
    b = bam_init1();
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        idx = sam_index_load(data[i]->fp, bam_name[i]);
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            std::shared_ptr<SAMRecord> read(new SAMRecord(b, data[i]->header));
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(read);
                }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(read);
                }
            }
        }
        hts_itr_destroy(iter);
    }

    bam_destroy1(b);
}

AlignmentContext ReadCache::getAlignmentContext() {
    getNextPos();
    std::vector<std::shared_ptr<SAMRecord>> tumor;
    std::vector<std::shared_ptr<SAMRecord>> normal;

    std::list<std::shared_ptr<SAMRecord>>::iterator iter = tumorReadsForAlignment.begin();
    while(iter != tumorReadsForAlignment.end()) {
        if((*iter)->getEndAfterFliter() < currentPose) {
            tumorReadsForAlignment.erase(iter++);
        } else {
            iter++;
        }
    }
    iter = normalReadsForAlignment.begin();
    while(iter != normalReadsForAlignment.end()) {
        if((*iter)->getEndAfterFliter() < currentPose) {
            normalReadsForAlignment.erase(iter++);
        } else {
            iter++;
        }
    }

    while(!tumorReads.empty() && tumorReads.front()->getStart() <= currentPose){
        tumorReadsForAlignment.emplace_back(tumorReads.front());
        tumorReadsForRegion.emplace_back(tumorReads.front());
        tumorReads.pop();
    }
    while(!normalReads.empty() && normalReads.front()->getStart() <= currentPose){
        normalReadsForAlignment.emplace_back(normalReads.front());
        normalReadsForRegion.emplace_back(normalReads.front());
        normalReads.pop();
    }
    SimpleInterval loc(data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceName(), currentPose, currentPose);
    for(std::shared_ptr<SAMRecord>& read : tumorReadsForAlignment) {
        if(!ReadUtils::isBaseInsideAdaptor(read, currentPose))
            tumor.emplace_back(read);
    }
    for(std::shared_ptr<SAMRecord>& read : normalReadsForAlignment) {
        if(!ReadUtils::isBaseInsideAdaptor(read, currentPose))
            normal.emplace_back(read);
    }
    return {tumor, normal, loc, tid, data[0]->header};
}

std::vector<std::shared_ptr<SAMRecord>> ReadCache::getReadsForRegion(AssemblyRegion & region) {
    std::vector<std::shared_ptr<SAMRecord>> ret;
    SimpleInterval loc = region.getExtendedSpan();
    std::list<std::shared_ptr<SAMRecord>>::iterator iter = tumorReadsForRegion.begin();
    while(iter != tumorReadsForRegion.end()) {
        SimpleInterval readLoc = (*iter)->getLoc();


        if((*iter)->getEndAfterFliter() < loc.getStart()) {
            tumorReadsForRegion.erase(iter++);
        } else {
            if(loc.overlaps(&readLoc)) {
                ret.emplace_back((*iter));
            }
            iter++;
        }

    }
    iter = normalReadsForRegion.begin();
    while(iter != normalReadsForRegion.end()) {
        SimpleInterval readLoc = (*iter)->getLoc();


        if((*iter)->getEndAfterFliter() < loc.getStart()) {
            normalReadsForRegion.erase(iter++);
        } else {
            if(loc.overlaps(&readLoc)) {
                ret.emplace_back((*iter));
            }
            iter++;
        }
    }
    return ret;
}
