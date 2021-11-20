//
// Created by 梦想家xixi on 2021/10/20.
//

#include "BaseVertex.h"
#include "Mutect2Utils.h"

int BaseVertex::hashCode(uint8_t *a, int length) {
    if(a == nullptr)
        return 0;
    int result = 1;
    for(int i = 0; i < length; i++) {
        result = 31 * result + a[i];
    }

    return result;
}

BaseVertex::BaseVertex(uint8_t * const sequence, const int length) : sequence(sequence), length(length){
    Mutect2Utils::validateArg(sequence != nullptr, "Sequence cannot be null");
    cashedHashCode = hashCode(sequence, length);
}

bool BaseVertex::isEmpty() const {
    return length == 0;
}

bool BaseVertex::operator==(const BaseVertex &other) const {
    if(other.cashedHashCode != cashedHashCode || other.length != length)
        return false;
    for(int i = 0; i < length; i++)
        if(sequence[i] != other.sequence[i])
            return false;
    return true;
}

bool BaseVertex::operator<(const BaseVertex &other) const {
    if(length > other.length)
        return false;
    if(length == other.length || cashedHashCode > other.cashedHashCode)
        return false;
    for(int i = 0; i < length; i++)
        if(sequence[i] > other.sequence[i])
            return false;
    return true;
}

std::ostream & operator<<(std::ostream &os, const BaseVertex &baseVertex) {
    os << "baseVertex : ";
    for(int i = 0; i < baseVertex.length; i++)
        os << baseVertex.sequence[i];
    os << '.' << std::endl;
    return os;
}

void BaseVertex::setAdditionalInfo(const std::string &info) {
    additionalInfo = info;
}

bool BaseVertex::hasAmbiguousSequence() {
    for(int i = 0; i < length; i++) {
        uint8_t tmp = sequence[i];
        if(tmp > 60)
            tmp -= 32;
        switch (tmp) {
            case 'A':
            case 'T':
            case 'G':
            case 'C':
                continue;
            default:
                return true;
        }
    }
    return false;
}

bool BaseVertex::seqEquals(BaseVertex *other) {
    if(length != other->getLength())
        return false;

    uint8_t * otherSeq = other->getSequence();
    for(int i = 0; i < length; i++){
        if(otherSeq[i] != sequence[i])
            return false;
    }
    return true;
}
