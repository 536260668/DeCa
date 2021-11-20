//
// Created by 梦想家xixi on 2021/11/15.
//

#include "SeqVertex.h"

long SeqVertex::hashCode() const {
    return (long)this;
}

bool SeqVertex::operator<(const SeqVertex &other) const {
    return hashCode() < other.hashCode();
}

SeqVertex *SeqVertex::withoutSuffix(uint8_t *suffix, int length) {
    int prefixSize = getLength() - length;
    int newLength = std::min(prefixSize, getLength());
    if(prefixSize > 0) {
        uint8_t * newSequence = new uint8_t[newLength];
        memcpy(newSequence, sequence, newLength);
        return new SeqVertex(newSequence, newLength);
    } else {
        return nullptr;
    }
}

SeqVertex *SeqVertex::withoutPrefixAndSuffix(uint8_t *prefix, int preLength, uint8_t *suffix, int sufLength) {
    int start = preLength;
    int length = getLength() - sufLength - preLength;
    int stop = start + length;
    int newLength = stop - start;
    if(length > 0) {
        uint8_t* newSequence = new uint8_t[newLength];
        memcpy(newSequence, sequence+start, newLength);
        return new SeqVertex(newSequence, newLength);
    } else {
        return nullptr;
    }
}
