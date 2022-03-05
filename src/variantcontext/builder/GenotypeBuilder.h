//
// Created by 梦想家xixi on 2021/11/30.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPEBUILDER_H
#define MUTECT2CPP_MASTER_GENOTYPEBUILDER_H

#include <utility>

#include "FastGenotype.h"
class GenotypeBuilder {
private:
    static std::vector<std::shared_ptr<Allele>> HAPLOID_NO_CALL;
    static std::vector<std::shared_ptr<Allele>> DIPLOID_NO_CALL;
    std::string sampleName;
    std::vector<std::shared_ptr<Allele>> alleles;
    bool isPhased = false;
    int GQ = -1;
    int DP = -1;
    int* AD = nullptr;
    int ADLength = 0;
    int* PL = nullptr;
    int PLLength = 0;
    std::map<std::string, void*> extendedAttributes;
    std::string filters;
    int initialAttributeMapSize = 5;
    static std::map<std::string, void*> NO_ATTRIBUTES;

public:
    static Genotype* create(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles);
    static Genotype* create(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles, const std::map<std::string, void*>& attributes);
    static Genotype* create(const std::string& sampleName, const std::vector<std::shared_ptr<Allele>>& alleles, double * gls, int length);
    GenotypeBuilder(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles) : sampleName(std::move(sampleName)), alleles(std::move(alleles)){}
    GenotypeBuilder attributes(const std::map<std::string, void*>& attributes);
    GenotypeBuilder buildPL(double * GLs, int length);
    Genotype* make();
};


#endif //MUTECT2CPP_MASTER_GENOTYPEBUILDER_H
