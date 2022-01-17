//
// Created by 梦想家xixi on 2021/12/1.
//

#include "AssemblyResultSet.h"

#include <utility>
#include "param/ParamUtils.h"

bool AssemblyResultSet::add(std::shared_ptr<Haplotype> &h, std::shared_ptr<AssemblyResult> &ar) {
    Mutect2Utils::validateArg(h.get(), "input haplotype cannot be null");
    Mutect2Utils::validateArg(ar.get(), "input assembly-result cannot be null");
    Mutect2Utils::validateArg(h->getGenomeLocation(), "the haplotype provided must have a genomic location");

    bool assemblyResultAdditionReturn = add(ar);

    if(haplotypes.find(h) != haplotypes.end()) {
        std::shared_ptr<AssemblyResult> previousAr = assemblyResultByHaplotype.at(h);
        if(previousAr == nullptr) {
            assemblyResultByHaplotype.insert(std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<AssemblyResult>>(h, ar));
            return true;
        } else if (previousAr != ar) {
            throw std::invalid_argument("there is already a different assembly result for the input haplotype");
        } else {
            return assemblyResultAdditionReturn;
        }
    } else {
        haplotypes.insert(h);
        assemblyResultByHaplotype.insert(std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<AssemblyResult>>(h, ar));
        updateReferenceHaplotype(h);
        if(h->getIsNonReference()) {
            variationPresent = true;
        }
        return true;
    }
}

bool AssemblyResultSet::add(std::shared_ptr<AssemblyResult> &ar) {
    Mutect2Utils::validateArg(ar.get(), "input assembly-result cannot be null");
    int kmerSize = ar->getKmerSize();
    if(assemblyResultByKmerSize.find(kmerSize) != assemblyResultByKmerSize.end()) {
        if(assemblyResultByKmerSize.at(kmerSize) != ar) {
            throw std::invalid_argument("a different assembly result with the same kmerSize was already added");
        }
        return false;
    } else {
        assemblyResultByKmerSize.insert(std::pair<int, std::shared_ptr<AssemblyResult>>(kmerSize, ar));
        kmerSizes.insert(kmerSize);
        return true;
    }
}

void AssemblyResultSet::updateReferenceHaplotype(std::shared_ptr<Haplotype> &newHaplotype) {
    if(!newHaplotype->getIsReference()) {
        return;
    }
    if(refHaplotype == nullptr) {
        refHaplotype = newHaplotype;
    } else {
        throw  std::invalid_argument("the assembly-result-set already have a reference haplotype that is different");
    }
}

void AssemblyResultSet::setRegionForGenotyping(AssemblyRegion & regionForGenotyping) {
    this->regionForGenotyping = &regionForGenotyping;
}

void AssemblyResultSet::setFullReferenceWithPadding(uint8_t *fullReferenceWithPadding, int length) {
    this->fullReferenceWithPadding = fullReferenceWithPadding;
    this->fullReferenceWithPaddingLength = length;
}

void AssemblyResultSet::setPaddedReferenceLoc(SimpleInterval *paddedReferenceLoc) {
    this->paddedReferenceLoc = paddedReferenceLoc;
}

bool AssemblyResultSet::add(std::shared_ptr<Haplotype> &h) {
    Mutect2Utils::validateArg(h.get(), "input haplotype can not be null");
    Mutect2Utils::validateArg(h->getGenomeLocation(), "haplotype genomeLocation cannot be null");
    if(haplotypes.find(h) != haplotypes.end()) {
        return false;
    }
    haplotypes.insert(h);
    updateReferenceHaplotype(h);
    return true;
}

std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & AssemblyResultSet::getVariationEvents(int maxMnpDistance) {
    ParamUtils::isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
    bool sameMnpDistance = lastMaxMnpDistanceUsed != -1 && maxMnpDistance == lastMaxMnpDistanceUsed;
    lastMaxMnpDistanceUsed = maxMnpDistance;
    bool flag = false;
    for(const std::shared_ptr<Haplotype>& haplotype : haplotypes) {
        if(haplotype->getIsReference() && haplotype->getEventMap()->empty()) {
            flag = true;
            break;
        }
    }
    if(variationEvents.empty() || !sameMnpDistance) {
        regenerateVariationEvents(maxMnpDistance);
    }
    return variationEvents;
}

void AssemblyResultSet::regenerateVariationEvents(int maxMnpDistance) {
    std::vector<std::shared_ptr<Haplotype>> haplotypeList = getHaplotypeList();
    EventMap::buildEventMapsForHaplotypes(haplotypeList, fullReferenceWithPadding, fullReferenceWithPaddingLength, paddedReferenceLoc, false, maxMnpDistance);
    variationEvents = EventMap::getAllVariantContexts(haplotypeList);
    lastMaxMnpDistanceUsed = maxMnpDistance;
    for(std::shared_ptr<Haplotype> haplotype : haplotypeList) {
        if(haplotype->getIsNonReference()) {
            variationPresent = true;
            break;
        }
    }
}

std::vector<std::shared_ptr<Haplotype>> AssemblyResultSet::getHaplotypeList() {
    return {haplotypes.begin(), haplotypes.end()};
}
