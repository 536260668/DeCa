//
// Created by 梦想家xixi on 2021/11/30.
//

#ifndef MUTECT2CPP_MASTER_FASTGENOTYPE_H
#define MUTECT2CPP_MASTER_FASTGENOTYPE_H

#include "Genotype.h"

class FastGenotype : public Genotype{
private:
    std::vector<std::shared_ptr<Allele>> alleles;
    bool _isPhased;
    int GQ;
    int DP;
    int* AD;
    int* PL;
    int ADLength;
    int PLLength;
    std::map<std::string, void*> extendedAttributes;

public:
    std::vector<std::shared_ptr<Allele>> & getAlleles() override {return alleles;}
    std::shared_ptr<Allele> getAllele(int i) override {return alleles.at(i);}
    bool isPhased() override {return _isPhased;}
    int getDP() override {return DP;}
    int* getAD(int& length) override;
    int getGQ() override {return GQ;}
    int* getPL(int & length) override;
    std::map<std::string, void*>& getExtendedAttributes() override {return extendedAttributes;}
    FastGenotype(std::string sampleName, std::vector<std::shared_ptr<Allele>> & alleles, bool isPhased, int GQ, int DP, int* AD, int ADLength, int* PL, int PLLength, const std::string& filters, std::map<std::string, void*> extendedAttributes);
};


#endif //MUTECT2CPP_MASTER_FASTGENOTYPE_H
