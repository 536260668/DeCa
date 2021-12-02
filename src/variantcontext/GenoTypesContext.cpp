//
// Created by lhh on 11/12/21.
//

#include "GenoTypesContext.h"

#include <utility>


GenoTypesContext::GenoTypesContext(int n) : maxPloidy(-1), immutable(false),  notToBeDirectlyAccessedGenotypes(std::vector<Genotype*>(n)){}

GenoTypesContext::GenoTypesContext(std::vector<Genotype *> genotypes) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(std::move(genotypes)){}
