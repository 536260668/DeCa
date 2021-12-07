//
// Created by lhh on 11/12/21.
//

#include "GenoTypesContext.h"

#include <utility>


GenoTypesContext::GenoTypesContext(int n) : maxPloidy(-1), immutable(false),  notToBeDirectlyAccessedGenotypes(new std::vector<Genotype*>(n)){}

GenoTypesContext::GenoTypesContext(std::vector<Genotype *> & genotypes) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(&genotypes), sampleNameToOffset(
        nullptr), sampleNamesInOrder(nullptr){}

GenoTypesContext GenoTypesContext::NO_GENOTYPES = GenoTypesContext(new std::vector<Genotype*>(), new std::map<std::string, int>(), new std::vector<std::string>()).setImmutable();

GenoTypesContext::GenoTypesContext(std::vector<Genotype *> *genotypes, std::map<std::string, int> *sampleNameToOffset,
                                   std::vector<std::string> *sampleNamesInOrder) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(genotypes), sampleNamesInOrder(sampleNamesInOrder), sampleNameToOffset(sampleNameToOffset){}

GenoTypesContext & GenoTypesContext::setImmutable() {
    immutable = true;
    return *this;
}

std::vector<Genotype *> *GenoTypesContext::getGenotypes() {
    return notToBeDirectlyAccessedGenotypes;
}

int GenoTypesContext::getSize() {
    return (int)getGenotypes()->size();
}

Genotype *GenoTypesContext::get(int i) {
    return getGenotypes()->at(i);
}
