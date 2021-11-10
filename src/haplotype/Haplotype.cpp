//
// Created by 梦想家xixi on 2021/11/9.
//

#include "Haplotype.h"
#include "read/AlignmentUtils.h"

Haplotype::Haplotype(uint8_t *bases, int length, bool isRef) : Allele(copyArray(bases, length), length, isRef){}

uint8_t *Haplotype::copyArray(uint8_t *base, int length) {
    uint8_t * res = new uint8_t[length];
    memcpy(res, base, length);
    return res;
}

Haplotype::Haplotype(uint8_t *bases, int length) : Allele(copyArray(bases, length), length, false){}

void Haplotype::setCigar(Cigar *cigar) {
    this->cigar = *AlignmentUtils::consolidateCigar(&cigar);
    Mutect2Utils::validateArg(this->cigar.getReadLength() == getLength(), "Read length is not equal to the read length of the cigar");
}

Haplotype::Haplotype(uint8_t *bases, bool isRef, int length, int alignmentStartHapwrtRef, Cigar *cigar) : Allele(copyArray(bases, length), length, false), alignmentStartHapwrtRef(alignmentStartHapwrtRef){
    setCigar(cigar);
}
