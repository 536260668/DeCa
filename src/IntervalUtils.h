//
// Created by 梦想家xixi on 2021/10/12.
//

#ifndef MUTECT2CPP_MASTER_INTERVALUTILS_H
#define MUTECT2CPP_MASTER_INTERVALUTILS_H

#include "SimpleInterval.h"

class SimpleInterval;

class IntervalUtils {
public:
    /**
     * Create a new interval, bounding start and stop by the start and end of contig
     *
     * This function will return null if start and stop cannot be adjusted in any reasonable way
     * to be on the contig.  For example, if start and stop are both past the end of the contig,
     * there's no way to fix this, and null will be returned.
     *
     * @param contig our contig
     * @param start our start as an arbitrary integer (may be negative, etc)
     * @param stop our stop as an arbitrary integer (may be negative, etc)
     * @param contigLength length of the contig
     * @return a valid interval over contig, or null if a meaningful interval cannot be created
     */
    static SimpleInterval* trimIntervalToContig(std::string contig, int start, int stop, int contigLength);

};

#endif //MUTECT2CPP_MASTER_INTERVALUTILS_H

