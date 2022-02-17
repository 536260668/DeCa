/**
 * The implementation of ReferenceCache class
 */
#include <iostream>
#include "ReferenceCache.h"
#include "assert.h"


ReferenceCache::ReferenceCache(char * refName, SAMFileHeader* header) : tid(0), header(header)
{
    fai = fai_load3_format(refName, NULL, NULL, FAI_CREATE, FAI_FASTA);
    start = 0;
    end = std::min(999999, header->getSequenceDictionary().getSequences()[tid].getSequenceLength());
    std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start) + '-' + std::to_string(end);
    hts_pos_t seq_len;
    bases = fai_fetch64(fai, region.c_str(), &seq_len);
}

ReferenceCache::~ReferenceCache()
{
    clear();
    fai_destroy(fai);
}

std::string ReferenceCache::getContig()
{
    return header->getSequenceDictionary().getSequences()[tid].getSequenceName();
}

void ReferenceCache::clear()
{
    if (this->bases && strlen(this->bases)) // strlen does not check for null pointer
        free(bases);
}




char ReferenceCache::getBase(hts_pos_t pos)
{
    while(pos > end) {
        advanceLoad();
    }
    return bases[pos - start];
}

void ReferenceCache::advanceLoad() {
    start = end + 1;
    end = std::max(start + 999999, static_cast<hts_pos_t>(header->getSequenceDictionary().getSequences()[tid].getSequenceLength()));
    std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start) + '-' + std::to_string(end);
    hts_pos_t seq_len;
    clear();
    bases = fai_fetch64(fai, region.c_str(), &seq_len);;
}

void ReferenceCache::setTid(int tid) {
    start = 0;
    end = std::max(99999, header->getSequenceDictionary().getSequences()[tid].getSequenceLength());
    std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start) + '-' + std::to_string(end);
    hts_pos_t seq_len;
    clear();
    bases = fai_fetch64(fai, region.c_str(), &seq_len);
}

std::shared_ptr<uint8_t> ReferenceCache::getSubsequenceAt(int tid, int start, int stop, int & length) {
    if(tid == this->tid && start >= this->start && stop <= this->end) {
        std::shared_ptr<uint8_t> ret(new uint8_t[stop-start+1]{0});
        std::copy(bases+start, bases+stop + 1, ret.get());
        length = stop - start + 1;
        return ret;
    }
    else {
        std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start) + '-' + std::to_string(stop);
        hts_pos_t seq_len;
        uint8_t * ret = reinterpret_cast<uint8_t*>(fai_fetch64(fai, region.c_str(), &seq_len));
        std::shared_ptr<uint8_t> toRet(new uint8_t[seq_len]);
        std::copy(ret, ret + seq_len, toRet.get());
        length = seq_len;
        return toRet;
    }
}
