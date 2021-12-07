//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
#define MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H

#include "Genotype.h"

class Genotype;

class GenoTypesContext {
protected:
    std::vector<std::string> * sampleNamesInOrder;
    std::map<std::string, int> * sampleNameToOffset;
    std::vector<Genotype*> * notToBeDirectlyAccessedGenotypes;
    std::vector<Genotype*> * getGenotypes();

private:
    int maxPloidy;
    bool immutable;
    std::vector<Genotype*> bases;

public:
    explicit GenoTypesContext(int n = 10);
    explicit GenoTypesContext(std::vector<Genotype*> & genotypes);
    GenoTypesContext(std::vector<Genotype*> * genotypes, std::map<std::string, int> * sampleNameToOffset, std::vector<std::string> * sampleNamesInOrder);
    static GenoTypesContext NO_GENOTYPES;
    GenoTypesContext & setImmutable();
    int getSize();
    Genotype* get(int i);
};


#endif //MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
