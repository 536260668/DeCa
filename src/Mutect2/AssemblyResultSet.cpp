//
// Created by 梦想家xixi on 2021/12/1.
//

#include "AssemblyResultSet.h"

bool AssemblyResultSet::add(Haplotype *h, AssemblyResult *ar) {
    Mutect2Utils::validateArg(h, "input haplotype cannot be null");
    Mutect2Utils::validateArg(ar, "input assembly-result cannot be null");
    Mutect2Utils::validateArg(h->getGenomeLocation(), "the haplotype provided must have a genomic location");

    bool assemblyResultAdditionReturn = add(ar);

    if(haplotypes.find(h) != haplotypes.end()) {
        AssemblyResult* previousAr = assemblyResultByHaplotype.at(h);
        if(previousAr == nullptr) {
            assemblyResultByHaplotype.insert(std::pair<Haplotype*, AssemblyResult*>(h, ar));
            return true;
        } else if (previousAr != ar) {
            throw std::invalid_argument("there is already a different assembly result for the input haplotype");
        } else {
            return assemblyResultAdditionReturn;
        }
    } else {
        haplotypes.insert(h);
        assemblyResultByHaplotype.insert(std::pair<Haplotype*, AssemblyResult*>(h, ar));
        updateReferenceHaplotype(h);
        if(h->getIsNonReference()) {
            variationPresent = true;
        }
        return true;
    }
}

bool AssemblyResultSet::add(AssemblyResult *ar) {
    Mutect2Utils::validateArg(ar, "input assembly-result cannot be null");
    int kmerSize = ar->getKmerSize();
    if(assemblyResultByKmerSize.find(kmerSize) != assemblyResultByKmerSize.end()) {
        if(assemblyResultByKmerSize.at(kmerSize) != ar) {
            throw std::invalid_argument("a different assembly result with the same kmerSize was already added");
        }
        return false;
    } else {
        assemblyResultByKmerSize.insert(std::pair<int, AssemblyResult*>(kmerSize, ar));
        kmerSizes.insert(kmerSize);
        return true;
    }
}

void AssemblyResultSet::updateReferenceHaplotype(Haplotype *newHaplotype) {
    if(!newHaplotype->getIsReference()) {
        return;
    }
    if(refHaplotype == nullptr) {
        refHaplotype = newHaplotype;
    } else {
        throw  std::invalid_argument("the assembly-result-set already have a reference haplotype that is different");
    }
}

void AssemblyResultSet::setRegionForGenotyping(AssemblyRegion *regionForGenotyping) {
    this->regionForGenotyping = regionForGenotyping;
}

void AssemblyResultSet::setFullReferenceWithPadding(uint8_t *fullReferenceWithPadding, int length) {
    this->fullReferenceWithPadding = fullReferenceWithPadding;
    this->fullReferenceWithPaddingLength = length;
}

void AssemblyResultSet::setPaddedReferenceLoc(SimpleInterval *paddedReferenceLoc) {
    this->paddedReferenceLoc = paddedReferenceLoc;
}

bool AssemblyResultSet::add(Haplotype *h) {
    Mutect2Utils::validateArg(h, "input haplotype can not be null");
    Mutect2Utils::validateArg(h->getGenomeLocation(), "haplotype genomeLocation cannot be null");
    if(haplotypes.find(h) != haplotypes.end()) {
        return false;
    }
    haplotypes.insert(h);
    updateReferenceHaplotype(h);
    return true;
}
