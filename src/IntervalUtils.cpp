//
// Created by 梦想家xixi on 2021/10/12.
//

#include "IntervalUtils.h"
#include "Mutect2Utils.h"

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

bool IntervalUtils::isAfter(Locatable &first, Locatable &second, SAMSequenceDictionary &dictionary) {
    int contigComparison = compareContigs(first, second, dictionary);
    return contigComparison == 1 || (contigComparison == 0 && first.getStart() > second.getEnd());
}

int IntervalUtils::compareContigs(Locatable &first, Locatable &second, SAMSequenceDictionary &dictionary) {
    std::string contig = first.getContig();
    int firstRefIndex = dictionary.getSequenceIndex(contig);
    contig = second.getContig();
    int secondRefIndex = dictionary.getSequenceIndex(contig);
    Mutect2Utils::validateArg(firstRefIndex != -1 && secondRefIndex != -1, "Can't do comparison because Locatables' contigs not found in sequence dictionary");
    return Mutect2Utils::Int_compare(firstRefIndex, secondRefIndex);
}
