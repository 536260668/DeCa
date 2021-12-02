//
// Created by 梦想家xixi on 2021/11/30.
//

#include "GenotypeBuilder.h"

#include <utility>

Genotype *GenotypeBuilder::create(std::string sampleName, std::vector<Allele *> alleles) {
    return GenotypeBuilder(std::move(sampleName), std::move(alleles)).make();
}

Genotype *GenotypeBuilder::make() {
    return new FastGenotype(sampleName, alleles, isPhased, GQ, DP, AD, ADLength, PL, PLLength, filters, extendedAttributes);
}

Genotype *GenotypeBuilder::create(std::string sampleName, std::vector<Allele *> alleles,
                                  const std::map<std::string, void *>& attributes) {
    return GenotypeBuilder(std::move(sampleName), std::move(alleles)).attributes(attributes).make();
}

GenotypeBuilder GenotypeBuilder::attributes(const std::map<std::string, void *>& attributes) {
    for(std::pair<std::string, void *> pairToAdd : attributes) {
        extendedAttributes.insert(pairToAdd);
    }
    return *this;
}

Genotype *GenotypeBuilder::create(const std::string& sampleName, const std::vector<Allele *>& alleles, double *gls, int length) {
    return GenotypeBuilder(sampleName, alleles).buildPL(gls, length).make();
}

GenotypeBuilder GenotypeBuilder::buildPL(double *GLs , int length) {
    PL = GenotypeLikelihoods::fromLog10Likelihoods(GLs, length)->getAsPLs();
    return *this;
}
