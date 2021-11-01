//
// Created by 梦想家xixi on 2021/10/19.
//

#include "Kmer.h"
#include "Mutect2Utils.h"
#include <iostream>

Kmer::Kmer(uint8_t *kmer, const int length)  : bases(kmer), start(0), length(length){
    Mutect2Utils::validateArg(start >= 0, "start must be >= 0");
    Mutect2Utils::validateArg(length >= 0, "length must be >= 0");
    this->hash = hashCode(bases, start, length);
}

Kmer::Kmer(uint8_t *kmer, const int start, const int length) : bases(kmer), start(start), length(length){
    Mutect2Utils::validateArg(start >= 0, "start must be >= 0");
    Mutect2Utils::validateArg(length >= 0, "length must be >= 0");
    this->hash = hashCode(bases, start, length);
}

Kmer::Kmer(const Kmer &kmer) : bases(kmer.bases), start(kmer.start), length(kmer.length) {
    Mutect2Utils::validateArg(start >= 0, "start must be >= 0");
    Mutect2Utils::validateArg(length >= 0, "length must be >= 0");
    this->hash = hashCode(bases, start, length);
}

int Kmer::hashCode(const uint8_t *bases, const int start, const int length) {
    if(length == 0) {
        return 0;
    }
    int h = 0;
    for(int i = start, stop = start + length; i < stop; i++) {
        h = 31 * h + bases[i];
    }
    return h;
}

Kmer Kmer::subKmer(const int newStart, const int newLength) {
    Kmer subkmer = Kmer(bases, start + newStart, newLength);
    return subkmer;
}

uint8_t *Kmer::getBases() const {
    uint8_t * res = new uint8_t[length];
    memcpy(res, bases+start, sizeof(uint8_t)*length);
    return res;
}

int Kmer::getDifferingPositions(Kmer other, int maxDistance, int *differingIndeces, uint8_t *differingBases) {
    Mutect2Utils::validateArg(differingIndeces != nullptr, "Null object is not allowed here.");
    Mutect2Utils::validateArg(differingBases != nullptr, "Null object is not allowed here.");
    Mutect2Utils::validateArg(maxDistance > 0, "maxDistance must be positive");
    int dist = 0;
    if(length == other.length) {
        uint8_t * f2 = other.getBases();
        for(int i = 0; i < length; i++) {
            if(bases[start + i] != f2[i]) {
                differingIndeces[dist] = i;
                differingBases[dist++] = f2[i];
                if(dist > maxDistance) {
                    return -1;
                }
            }
        }
        delete[] f2;
    }
    return dist;
}

bool Kmer::operator<(const Kmer &other) const {
    if( this->length > other.length)
        return false;
    if(this->hash > other.hash)
        return false;
    for(int i = 0; i < length; i++)
        if(this->bases[this->start+i] > other.bases[other.start+i])
            return false;
    return true;
}

bool Kmer::operator==(const Kmer &other) const {
    if(this->hash != other.hash || this->length != other.length)
        return false;
    for(int i = 0; i < length; i++)
        if(this->bases[this->start+i] != other.bases[other.start+i])
            return false;
    return true;
}
