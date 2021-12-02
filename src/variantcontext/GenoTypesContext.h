//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
#define MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H

#include "Genotype.h"

class Genotype;

class GenoTypesContext {
protected:
    std::vector<std::string> sampleNamesInOrder;
    std::map<std::string, int> sampleNameToOffset;
    std::vector<Genotype*> notToBeDirectlyAccessedGenotypes;

private:
    int maxPloidy;
    bool immutable;
    std::vector<Genotype*> bases;

public:
    explicit GenoTypesContext(int n = 10);
    explicit GenoTypesContext(std::vector<Genotype*> genotypes);
};


#endif //MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
