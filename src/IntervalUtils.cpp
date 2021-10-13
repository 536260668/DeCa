//
// Created by 梦想家xixi on 2021/10/12.
//

#include "IntervalUtils.h"

SimpleInterval* IntervalUtils::trimIntervalToContig(const std::string contig, const int start, const int stop, const int contigLength) {
    if(contig.empty())
        throw std::invalid_argument("Null object is not allowed here.");
    if(contigLength < 1)
        throw std::invalid_argument("ContigLength should be at least 1.");

    const int boundedStart = std::max(1, start);
    const int boundedStop = std::min(contigLength, stop);

    if ( boundedStart > contigLength || boundedStop < 1 ){
        // there's no meaningful way to create this interval, as the start and stop are off the contig
        return nullptr;
    } else {
        return new SimpleInterval(contig, boundedStart, boundedStop);
    }
}