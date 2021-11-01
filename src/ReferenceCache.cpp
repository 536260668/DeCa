/**
 * The implementation of ReferenceCache class
 */
#include <iostream>
#include "ReferenceCache.h"
#include "assert.h"


ReferenceCache::ReferenceCache(char * refName)
{
    this->fai = fai_load3_format(refName, NULL, NULL, FAI_CREATE, FAI_FASTA);
    assert(fai != NULL);
    memset(this->chroName, '\0', sizeof(this->chroName));

    this->bases = NULL;
    this->start = -1;
    this->end = -1;
}

ReferenceCache::~ReferenceCache()
{
    clear();
    fai_destroy(fai);
}

char * ReferenceCache::getContig()
{
    return this->chroName;
}

void ReferenceCache::clear()
{
    if (this->bases && strlen(this->bases)) // strlen does not check for null pointer
        free(bases);
}

void ReferenceCache::fill(const char *chrmName, hts_pos_t start)
{
    char region[20];
    sprintf(region, "%s:%ld-", chrmName, start+1);  // region requires 1-based?
    //printf("start: %ld\n", start);
    hts_pos_t seq_len;
    bases = fai_fetch64(fai, region, &seq_len);
    assert(seq_len > 0);

    strcpy(this->chroName, chrmName);
    this->start = start;
    this->end = start + seq_len - 1;
}

// get sequence using with copying
const char * ReferenceCache::getSequence(hts_pos_t start, hts_pos_t end)
{
    assert(bases);  // bases shouldn't be NULL, otherwise strlen will throw an error
    if (!strlen(bases) || start > this->end)
        return NULL;

    int readLength = end - start + 1;
    char * refBases = new char[readLength + 1]; // rlen + '\0'

    if (start < this->start)    // if the read is not ordered
    {
        hts_pos_t fai_ref_len;
        //printf("read is not ordered\n");
	    return faidx_fetch_seq64(fai, this->chroName, start, end, &fai_ref_len);
    }
    hts_pos_t index = start - this->start;

    memcpy(refBases, bases + index, readLength);
    refBases[readLength] = '\0';

    //printf("%s\n", refBases);

    return refBases;
}

char ReferenceCache::getBase(hts_pos_t pos)
{
    assert(bases);  // bases shouldn't be NULL, otherwise strlen will throw an error
    if (!strlen(bases) || pos > this->end)
        return (char)NULL;

    if (pos < this->start)    // if the read is not ordered
    {
        hts_pos_t fai_ref_len;
        //printf("read is not ordered\n");
        char * base = faidx_fetch_seq64(fai, this->chroName, pos, pos, &fai_ref_len);
        char temp = *base;
        free(base);
        return temp;
    }

    hts_pos_t index = pos - this->start;
    return *(bases + index);
}
