//
// Created by 梦想家xixi on 2021/10/20.
//

#include "MultiDeBruijnVertex.h"

MultiDeBruijnVertex::MultiDeBruijnVertex(uint8_t *sequence, int length, bool mergeIdenticalNodes) : BaseVertex(sequence, length), mergeIdenticalNodes(mergeIdenticalNodes){
    hashCode = mergeIdenticalNodes ? BaseVertex::getHashCode() : (long) this;
}

bool MultiDeBruijnVertex::operator==(const MultiDeBruijnVertex &other) const {
    if(this->getLength() != other.getLength() || this->getHashCode() != other.getHashCode())
        return false;
    for(int i = 0; i < this->getLength(); i++)
        if(sequence[i] != other.sequence[i])
            return false;
    return true;
}

bool MultiDeBruijnVertex::operator<(const MultiDeBruijnVertex &other) const {
    if(this->getLength() > other.getLength())
        return false;
    if(this->getLength() == other.getLength() || this->getHashCode() > other.getHashCode())
        return false;
    for(int i = 0; i < this->getLength(); i++)
        if(sequence[i] > other.sequence[i])
            return false;
    return true;
}

uint8_t *MultiDeBruijnVertex::getAdditionalSequence(bool source) {
    return source ? BaseVertex::getAdditionalSequence(source) : getSuffixAsArray();
}

uint8_t *MultiDeBruijnVertex::getSuffixAsArray() const {
    uint8_t * res = new uint8_t[1];
    res[0] = getSuffix();
    return res;
}

MultiDeBruijnVertex::MultiDeBruijnVertex(uint8_t *sequence, int length) : BaseVertex(sequence, length) , mergeIdenticalNodes(false){
    hashCode = mergeIdenticalNodes ? BaseVertex::getHashCode() : (long) this;
}

int MultiDeBruijnVertex::getAdditionalLength(bool source) {
    return source ? getLength() : 1;
}

int MultiDeBruijnVertex::getAdditionalSequenceLength(bool isSource) {
    return isSource ? getLength() : 1;
}
