//
// Created by lhh on 11/11/21.
//

#ifndef MUTECT2CPP_MASTER_VARIANTCONTEXT_H
#define MUTECT2CPP_MASTER_VARIANTCONTEXT_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "htslib/sam.h"
#include "haplotype/Allele.h"
#include "VCFConstants.h"
#include "GenoTypesContext.h"
#include "CommonInfo.h"

//TODO: add GenoTypesContext class 2021.11.11
class VariantContext {
private:
    std::string contig;
    hts_pos_t start;
    hts_pos_t stop;
    std::string ID;

    CommonInfo commonInfo;

    /** A set of the alleles segregating in this context */
    std::vector<Allele> * alleles;


    std::vector<Allele> * makeAlleles(std::vector<Allele> & alleles);

public:
    /**
     * the actual constructor.  Private access only
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param genotypes       genotypes map
     * @param log10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @param validationToPerform     set of validation steps to take
     */
protected:
    VariantContext(std::string source,
                         std::string ID,
                         std::string contig,
                         long start,
                         long stop,
                         std::vector<Allele> alleles,
                         GenoTypesContext genotypes,
                         double log10PError,
                         std::set<std::string> filters     //TODO: finish the parameter list
                          );


};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXT_H
