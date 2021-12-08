//
// Created by 梦想家xixi on 2021/12/6.
//

#ifndef MUTECT2CPP_MASTER_VARIANTCONTEXTBUILDER_H
#define MUTECT2CPP_MASTER_VARIANTCONTEXTBUILDER_H

#include <string>
#include <vector>
#include <set>
#include "Allele.h"
#include "variantcontext/GenoTypesContext.h"

class VariantContextBuilder {
private:
    bool fullyDecoded = false;
    std::string source;
    std::string contig;
    long start = -1;
    long stop = -1;
    std::vector<Allele*> * alleles;
    std::string ID = ".";
    GenoTypesContext* genotypes;
    double log10PError;
    std::set<std::string> * filters;
    std::map<std::string, void*> * attribute;
    bool attributesCanBeModified;
    std::set<Validation> toValidate;

public:
    VariantContextBuilder(std::string & source, std::string & contig, long start, long stop, std::vector<Allele*> * alleles);

    VariantContextBuilder(VariantContext * parent);

    VariantContext* make(bool leaveModifyableAsIs);

    VariantContext* make();

    void setStop(long stop);

    void setAlleles(std::vector<Allele*> * alleles);
};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXTBUILDER_H
