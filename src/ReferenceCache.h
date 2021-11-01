/**
 * A helper class to cache the reference sequence
 */

#ifndef REFERENCE_CACHE_H
#define REFERENCE_CACHE_H

#include "htslib/faidx.h"

class ReferenceCache
{
private:
    char * bases;
    char chroName[32];   // name of cached chromosome
    hts_pos_t start;
    hts_pos_t end;

    faidx_t * fai;
public:
    ReferenceCache(char * refName);

    ~ReferenceCache();

    // get the name of cached chromosome
    char * getContig();

    // delete all the elements in the cache
    void clear();

    void fill(const char *chrmName, hts_pos_t start);

    /**
     * Get sequence of a specific interval, without copying to another array
     */
    const char * getSequence(hts_pos_t start, hts_pos_t end);

    /**
     * Get a single base from the reference cache
     */
    char getBase(hts_pos_t pos);
};


#endif
