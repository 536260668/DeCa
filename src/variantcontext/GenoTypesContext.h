//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
#define MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
#include <unordered_map>
#include "Genotype.h"

class Genotype;

class GenoTypesContext {
protected:
    std::vector<std::string> * sampleNamesInOrder;
    std::unordered_map<std::string, int> * sampleNameToOffset;
    std::vector<Genotype*> * notToBeDirectlyAccessedGenotypes;


    void ensureSampleNameMap();
    void invalidateSampleOrdering();

private:
    int maxPloidy;
    bool immutable = false;
    std::vector<Genotype*> bases;

public:
    static GenoTypesContext NO_GENOTYPES;

    explicit GenoTypesContext(int n = 10);
    explicit GenoTypesContext(std::vector<Genotype*> & genotypes);
    GenoTypesContext(std::vector<Genotype*> * genotypes, std::unordered_map<std::string, int> * sampleNameToOffset, std::vector<std::string> * sampleNamesInOrder);
    ~GenoTypesContext();

    std::vector<Genotype*> * getGenotypes();
    GenoTypesContext & setImmutable();
    int getSize();
    Genotype* get(int i);
    bool containsSample(std::string sample);
    bool add(Genotype* genotype);
    bool isEmpty();
};


#endif //MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
