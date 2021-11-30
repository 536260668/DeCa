//
// Created by lhh on 11/11/21.
//

#include "VariantContext.h"

VariantContext::VariantContext(std::string source, std::string ID, std::string contig, long start, long stop,
                               std::vector<Allele> alleles, GenoTypesContext genotypes, double log10PError,
                               std::set<std::string> filters) : contig(contig), start(start), stop(stop), commonInfo(source, log10PError, filters)
{
    if(ID.empty() || std::equal(ID.begin(), ID.end(), ""))
        throw "ID field cannot be the null or the empty string";

    this->ID = ID == VCFConstants::EMPTY_ID_FIELD ? VCFConstants::EMPTY_ID_FIELD : ID;

    this->alleles = makeAlleles(alleles);

    // TODO: finish this method 2021.11.13
}

//std::vector<Allele> *VariantContext::makeAlleles(std::vector<Allele> &alleles)
//{
//    std::vector<Allele> * alleleList = new std::vector<Allele>(alleles.size());
//    bool sawRef = false;
//    for(Allele a : alleles)
//    {
//        for(int i=0, alleleListSize = alleleList->size(); i<alleleListSize; i++)
//        {
//            //---TODO: how to judge whether two allele object are equal ?
//        }
//
//        // deal with the case where the first allele isn't the reference
//        if(a.getIsReference())
//        {
//            if(sawRef)
//                throw "Alleles for a VariantContext must contain at most one reference allele: ";
//            alleleList->insert(alleleList->begin(), a);
//            sawRef = true;
//        } else {
//            alleleList->push_back(a);
//        }
//    }
//
//    if(alleleList->empty())
//        throw "Cannot create a VariantContext with an empty allele list";
//
//    if(alleleList->at(0).getIsNonReference())
//        throw "Alleles for a VariantContext must contain at least one reference allele: ";
//
//    return alleleList;
//}