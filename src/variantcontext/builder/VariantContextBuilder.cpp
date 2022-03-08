//
// Created by 梦想家xixi on 2021/12/6.
//

#include "VariantContextBuilder.h"

VariantContextBuilder::VariantContextBuilder(std::string & source, std::string & contig, long start, long stop,
                                             const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> &  alleles) : genotypes(&GenoTypesContext::NO_GENOTYPES), log10PError(1.0),
                                             attributesCanBeModified(false), source(source), contig(contig), start(start), stop(stop), alleles(alleles), filters(
                nullptr), attribute(nullptr){
    toValidate.insert(ALLELES);
}

std::shared_ptr<VariantContext> VariantContextBuilder::make(bool leaveModifyableAsIs) {
    if(!leaveModifyableAsIs) {
        attributesCanBeModified = false;
    }

    std::shared_ptr<VariantContext> ret = std::make_shared<VariantContext>(source, ID, contig, start, stop, alleles, genotypes, log10PError, filters, attribute, fullyDecoded,toValidate);
    return ret;
}

std::shared_ptr<VariantContext> VariantContextBuilder::make() {
    return make(false);
}

VariantContextBuilder::VariantContextBuilder(std::shared_ptr<VariantContext> &parent) : alleles(std::make_shared<std::vector<std::shared_ptr<Allele>>>(parent->getAlleles()) ), attribute(&parent->getAttributes()),attributesCanBeModified(false),
contig(parent->getContig()), filters(parent->getFiltersMaybeNull()), genotypes(parent->getGenotypes()), ID(parent->getID()), log10PError(parent->getLog10PError()),
source(parent->getSource()), start(parent->getStart()), stop(parent->getEnd()), fullyDecoded(parent->isFullyDecoded()){}

void VariantContextBuilder::setStop(long stop) {
    this->stop = stop;
}

void VariantContextBuilder::setAlleles(const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> &  alleles) {
    this->alleles = alleles;
    toValidate.insert(ALLELES);
}

VariantContextBuilder::~VariantContextBuilder() {
    alleles->clear();
}
