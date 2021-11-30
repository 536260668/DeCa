//
// Created by 梦想家xixi on 2021/11/30.
//

#ifndef MUTECT2CPP_MASTER_FASTGENOTYPE_H
#define MUTECT2CPP_MASTER_FASTGENOTYPE_H

#include "Genotype.h"

class FastGenotype : public Genotype{
private:
    std::vector<Allele*> alleles;
    bool _isPhased;
    int GQ;
    int DP;
    int* AD;
    int* PL;
    int ADLength;
    int PLLength;
    std::map<std::string, void*> extendedAttributes;

public:
    std::vector<Allele*> & getAlleles() override {return alleles;}
    Allele* getAllele(int i) override {return alleles.at(i);}
    bool isPhased() override {return _isPhased;}

};


#endif //MUTECT2CPP_MASTER_FASTGENOTYPE_H
