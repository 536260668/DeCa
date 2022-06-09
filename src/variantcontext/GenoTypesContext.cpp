//
// Created by lhh on 11/12/21.
//

#include "GenoTypesContext.h"

#include <utility>
#include <unordered_map>
#include <cassert>


GenoTypesContext::GenoTypesContext(int n) : maxPloidy(-1), immutable(false),  notToBeDirectlyAccessedGenotypes(new std::vector<Genotype*>(n)), sampleNameToOffset(
        nullptr), sampleNamesInOrder(nullptr){}

GenoTypesContext::GenoTypesContext(std::vector<Genotype *> & genotypes) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(&genotypes), sampleNameToOffset(
        nullptr), sampleNamesInOrder(nullptr){}

        GenoTypesContext GenoTypesContext::NO_GENOTYPES = GenoTypesContext(nullptr, nullptr, nullptr).setImmutable();

GenoTypesContext::GenoTypesContext(std::vector<Genotype *> *genotypes, std::unordered_map<std::string, int> *sampleNameToOffset,
                                   std::vector<std::string> *sampleNamesInOrder) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(genotypes), sampleNamesInOrder(sampleNamesInOrder), sampleNameToOffset(sampleNameToOffset){}

GenoTypesContext::~GenoTypesContext() noexcept {
    if(notToBeDirectlyAccessedGenotypes)
    {
        for(Genotype* genotype : *notToBeDirectlyAccessedGenotypes)
        {
            delete genotype;
        }
    }

    delete notToBeDirectlyAccessedGenotypes;
    delete sampleNameToOffset;
    delete sampleNamesInOrder;
}

GenoTypesContext & GenoTypesContext::setImmutable() {
    immutable = true;
    return *this;
}

std::vector<Genotype *> *GenoTypesContext::getGenotypes() {
    return notToBeDirectlyAccessedGenotypes;
}

int GenoTypesContext::getSize() {
    if(!notToBeDirectlyAccessedGenotypes)
        return 0;
    else
        return (int)notToBeDirectlyAccessedGenotypes->size();
}

Genotype *GenoTypesContext::get(int i) {
    return getGenotypes()->at(i);
}

void GenoTypesContext::ensureSampleNameMap() {
    if(sampleNameToOffset == nullptr)
    {
        sampleNameToOffset = new std::unordered_map<std::string, int>(getSize());

        for(int i=0; i<getSize(); i++)
        {
            sampleNameToOffset->insert({notToBeDirectlyAccessedGenotypes->operator[](i)->getSampleName(), i});
        }
    }
}

bool GenoTypesContext::containsSample(std::string sample) {
    ensureSampleNameMap();
    return sampleNameToOffset->find(sample) != sampleNameToOffset->end();
}

bool GenoTypesContext::add(Genotype *genotype) {
    assert(!this->immutable);
    invalidateSampleOrdering();
    if(sampleNameToOffset != nullptr)
    {
        sampleNameToOffset->insert({genotype->getSampleName(), getSize()});
    }
    return getGenotypes()->emplace_back(genotype);
}

void GenoTypesContext::invalidateSampleOrdering() {
    sampleNamesInOrder = nullptr;
}

bool GenoTypesContext::isEmpty() {
    return notToBeDirectlyAccessedGenotypes == nullptr || notToBeDirectlyAccessedGenotypes->empty();
}