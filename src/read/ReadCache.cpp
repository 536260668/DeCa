//
// Created by 梦想家xixi on 2022/1/11.
//

#include "ReadCache.h"
#include "iostream"
#include "ReadUtils.h"


ReadCache::ReadCache(aux_t **data, std::vector<char*> & bam_name, std::shared_ptr<ReferenceCache> & cache) : data(data), tid(0), bam_name(bam_name),
                                                                    readTransformer(cache, data[0]->header, 5){
    bam1_t * b;
    b = bam_init1();
    std::string& region = data[0]->header->getSequenceDictionary().getSequences()[0].getSequenceName();
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        int count = 0;
        idx = sam_index_load(data[i]->fp, bam_name[i]);
        if(idx == 0)
            throw std::invalid_argument("random alignment retrieval only works for indexed BAM or CRAM files.");
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            std::shared_ptr<SAMRecord> read(new SAMRecord(b, data[i]->header));
            read = readTransformer.apply(read);
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(getpileRead(read));
                }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(getpileRead(read));
                }
                count++;
            }
            if(count > 500)
                break;
        }
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
    }
    bam_destroy1(b);

    start = end = currentPose = 0;
}

ReadCache::ReadCache(aux_t **data, std::vector<char *> &bam_name, int tid, const std::string& region, std::shared_ptr<ReferenceCache> & cache) : tid(tid), data(data), bam_name(bam_name), readTransformer(cache, data[0]->header, 5){
    bam1_t * b;
    b = bam_init1();
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        idx = sam_index_load(data[i]->fp, bam_name[i]);
        if(idx == 0)
            throw std::invalid_argument("random alignment retrieval only works for indexed BAM or CRAM files.");
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            std::shared_ptr<SAMRecord> read = std::make_shared<SAMRecord>(b, data[i]->header);
            read = readTransformer.apply(read);
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(getpileRead(read));
                    }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(getpileRead(read));
                }
            }
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);
    hts_idx_destroy(idx);
    unsigned i = region.find_last_of(':') + 1;
    start = end = 0;
    while(region[i] >= '0' && region[i] <= '9') {
        start = start * 10 + region[i] - '0';
        i++;
    }
    currentPose = start - 1;
    i = region.find_last_of('-') + 1;
    while(region[i] >= '0' && region[i] <= '9') {
        end = end * 10 + region[i] - '0';
        i++;
    }
}

ReadCache::~ReadCache() {
    for(pileRead* read : normalReadsForAlignment) {
        delete read;
    }
    for(pileRead* read : tumorReadsForAlignment) {
        delete read;
    }

    for(pileRead* read : normalCache) {
        delete read;
    }

    for(pileRead* read : tumorCache) {
        delete read;
    }
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
            std::shared_ptr<SAMRecord> read = std::make_shared<SAMRecord>(b, data[i]->header);
            if(ReadFilter::test(read, data[i]->header)) {
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(getpileRead(read));
                }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(getpileRead(read));
                }
            }
        }
        hts_itr_destroy(iter);
    }

    bam_destroy1(b);
}

AlignmentContext ReadCache::getAlignmentContext() {
    getNextPos();

    std::list<pileRead*>::iterator iter = tumorReadsForAlignment.begin();
    while(iter != tumorReadsForAlignment.end() && (*iter)->activateStop < currentPose) {
        delete *iter;
        tumorReadsForAlignment.erase(iter++);
    }
    iter = normalReadsForAlignment.begin();
    while(iter != normalReadsForAlignment.end() && (*iter)->activateStop < currentPose) {
        delete *iter;
        normalReadsForAlignment.erase(iter++);
    }

    while(!tumorReads.empty() && (*tumorReads.front()).read->getStart() <= currentPose){
        if((*tumorReads.front()).activateStart <= currentPose) {
            if((*tumorReads.front()).activateStop >= currentPose) {
                InsertPileToAlignment(tumorReads.front(), tumorReadsForAlignment);
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                tumorReads.pop();
            } else {
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                delete tumorReads.front();
                tumorReads.pop();
            }
        } else {
            if(InsertPileToCache(tumorReads.front(), tumorCache)) {
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                tumorReads.pop();
            } else {
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                delete tumorReads.front();
                tumorReads.pop();
            }
        }

    }
    while(!normalReads.empty() && (*normalReads.front()).read->getStart() <= currentPose){
        if((*normalReads.front()).activateStart <= currentPose) {
            if((*normalReads.front()).activateStop >= currentPose) {
                InsertPileToAlignment(normalReads.front(), normalReadsForAlignment);
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                normalReads.pop();
            } else {
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                delete normalReads.front();
                normalReads.pop();
            }
        } else {
            if(InsertPileToCache(normalReads.front(), normalCache)) {
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                normalReads.pop();
            } else {
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                delete normalReads.front();
                normalReads.pop();
            }
        }

    }
    iter = tumorCache.begin();
    while (iter != tumorCache.end() && (*iter)->activateStart <= currentPose) {
        InsertPileToAlignment(*iter, tumorReadsForAlignment);
        tumorCache.erase(iter++);
    }
    iter = normalCache.begin();
    while (iter != normalCache.end() && (*iter)->activateStart <= currentPose) {
        InsertPileToAlignment(*iter, normalReadsForAlignment);
        normalCache.erase(iter++);
    }
    SimpleInterval loc(data[0]->header->getSequenceDictionary().getSequences()[tid].getSequenceName(), currentPose, currentPose);
    return {tumorReadsForAlignment, normalReadsForAlignment, loc, tid, data[0]->header};
}

std::vector<std::shared_ptr<SAMRecord>> ReadCache::getReadsForRegion(AssemblyRegion & region) {
    std::vector<std::shared_ptr<SAMRecord>> ret;
    std::shared_ptr<SimpleInterval> loc = region.getExtendedSpan();
    std::list<std::shared_ptr<SAMRecord>>::iterator iter = tumorReadsForRegion.begin();
    while(iter != tumorReadsForRegion.end()) {
        std::shared_ptr<SimpleInterval> readLoc = (*iter)->getLoc();


        if((*iter)->getEndAfterFliter() < loc->getStart()) {
            tumorReadsForRegion.erase(iter++);
        } else {
            if(loc->overlaps(readLoc)) {
                ret.emplace_back((*iter));
            }
            iter++;
        }

    }
    iter = normalReadsForRegion.begin();
    while(iter != normalReadsForRegion.end()) {
        std::shared_ptr<SimpleInterval> readLoc = (*iter)->getLoc();


        if((*iter)->getEndAfterFliter() < loc->getStart()) {
            normalReadsForRegion.erase(iter++);
        } else {
            if(loc->overlaps(readLoc)) {
                ret.emplace_back((*iter));
            }
            iter++;
        }
    }
    return ret;
}

pileRead *ReadCache::getpileRead(const std::shared_ptr<SAMRecord> &read) {
    int adaptorBoundary = read->getAdaptorBoundary();
    if (adaptorBoundary == INT32_MIN || read->getFragmentLength() > 100)
        return new pileRead{read, read->getStart(), read->getEnd()};
    int start = 0;
    int end = 0;
    if(read->isReverseStrand()) {
        start = std::max(adaptorBoundary+1, read->getStart());
        end = read->getEnd();

    } else {
        start = read->getStart();
        end = std::min(adaptorBoundary-1, read->getEnd());
    }
    return new pileRead{read, start, end};
}

void ReadCache::InsertPileToAlignment(pileRead* pile, std::list<pileRead*> & toAdd) {
    std::list<pileRead*>::iterator iter = toAdd.end();
    iter--;
    while(iter != toAdd.begin())
    {
        if((*iter)->activateStop < pile->activateStop)
            break;
        iter--;
    }
    if(!toAdd.empty() && (*iter)->activateStop < pile->activateStop)
        toAdd.insert(++iter, pile);
    else {
        toAdd.emplace_front(pile);
    }
}

bool ReadCache::InsertPileToCache(pileRead *pile, std::list<pileRead *> & toAdd) {
    if(pile->activateStart > pile->activateStop)
        return false;
    std::list<pileRead*>::iterator iter = toAdd.end();
    iter--;
    while(iter != toAdd.begin())
    {
        if((*iter)->activateStart < pile->activateStart)
            break;
        iter--;
    }
    if(!toAdd.empty() && (*iter)->activateStart < pile->activateStart)
        toAdd.insert(++iter, pile);
    else {
        toAdd.emplace_front(pile);
    }
    return true;
}
