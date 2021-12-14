//
// Created by 梦想家xixi on 2021/11/9.
//

#include <cstring>
#include "Haplotype.h"
#include "read/AlignmentUtils.h"
#include "SimpleInterval.h"

Haplotype::Haplotype(uint8_t *bases, int length, bool isRef) : Allele(copyArray(bases, length), length, isRef){}

uint8_t *Haplotype::copyArray(uint8_t *base, int length) {
    uint8_t * res = new uint8_t[length];
    memcpy(res, base, length);
    return res;
}

Haplotype::Haplotype(uint8_t *bases, int length) : Allele(copyArray(bases, length), length, false){}

void Haplotype::setCigar(Cigar *cigar) {
    this->cigar = AlignmentUtils::consolidateCigar(cigar);
    Mutect2Utils::validateArg(this->cigar->getReadLength() == getLength(), "Read length is not equal to the read length of the cigar");
}

Haplotype::Haplotype(uint8_t *bases, bool isRef, int length, int alignmentStartHapwrtRef, Cigar *cigar) : Allele(copyArray(bases, length), length, false), alignmentStartHapwrtRef(alignmentStartHapwrtRef){
    setCigar(cigar);
}

Haplotype::Haplotype(uint8_t *bases, int length, Locatable *loc) : Allele(copyArray(bases, length), length, false), genomeLocation(loc){}

Haplotype *Haplotype::trim(Locatable* loc) {
    Mutect2Utils::validateArg(loc != nullptr, "Loc cannot be null");
    Mutect2Utils::validateArg(genomeLocation != nullptr, "Cannot trim a Haplotype without containing GenomeLoc");
    SimpleInterval interval = SimpleInterval(genomeLocation);
    Mutect2Utils::validateArg(interval.contains(loc), "Can only trim a Haplotype to a containing span.");
    Mutect2Utils::validateArg(getCigar() != nullptr, "Cannot trim haplotype without a cigar");

    int newStart = loc->getStart() - this->genomeLocation->getStart();
    int newStop = newStart + loc->getEnd() - loc->getStart();
    uint8_t * newBases = AlignmentUtils::getBasesCoveringRefInterval(newStart, newStop, getBases(), getBasesLength(), 0, getCigar());
    Cigar* newCigar = AlignmentUtils::trimCigarByReference(getCigar(), newStart, newStop);

    if(newBases == nullptr || AlignmentUtils::startsOrEndsWithInsertionOrDeletion(newCigar))
        return nullptr;

    Haplotype* ret = new Haplotype(newBases, getIsReference());
    ret->setCigar(newCigar);
    ret->setGenomeLocation(loc);
    ret->setScore(score);
    ret->setAlignmentStartHapwrtRef(newStart + getAlignmentStartHapwrtRef());
    return ret;
}

Cigar *Haplotype::getCigar() {
    return cigar;
}

void Haplotype::setGenomeLocation(Locatable *genomeLocation) {
    this->genomeLocation = genomeLocation;
}

void Haplotype::setScore(double score) {
    this->score = score;
}

int Haplotype::getAlignmentStartHapwrtRef() const{
    return alignmentStartHapwrtRef;
}

EventMap *Haplotype::getEventMap() {
    return eventMap;
}

void Haplotype::setEventMap(EventMap *eventMap) {
    this->eventMap = eventMap;
}

bool Haplotype::operator<(const Haplotype &other) const {
    if(this->getLength() < other.getLength())
        return true;
    else if(this->getLength() > other.getLength())
        return false;
    else {
        uint8_t * bases = this->getBases();
        uint8_t * otherBases = other.getBases();
        for(int i = 0; i < this->getLength(); i++) {
            if(bases[i] < otherBases[i])
                return true;
            else if (bases[i] > otherBases[i])
                return false;
            else {
                continue;
            }
        }
        return false;
    }
}


