//
// Created by lhh on 11/11/21.
//

#include <stdexcept>
#include <cassert>
#include "VariantContext.h"

VariantContext::VariantContext(std::string &source,
                               std::string &ID,
                               std::string &contig,
                               long start,
                               long stop,
                               const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> & alleles,
                               GenoTypesContext* genotypes,
                               double log10PError,
                               std::set<std::string>* filters, std::map<std::string, void*> *attributes,
                               bool fullyDecoded,
                               const std::set<Validation>&  validationToPerform) : contig(contig), start(start), stop(stop), commonInfo(
        CommonInfo(source, log10PError, filters))
{
    type = VariantContext_NULL;
    if(ID.empty() || std::equal(ID.begin(), ID.end(), ""))
        throw std::invalid_argument("ID field cannot be the null or the empty string");

    this->ID = ID == VCFConstants::EMPTY_ID_FIELD ? VCFConstants::EMPTY_ID_FIELD : ID;

    this->alleles = std::move(makeAlleles(*alleles));

    // TODO: finish this method 2021.11.13

    if(genotypes != nullptr && genotypes != &GenoTypesContext::NO_GENOTYPES) {
        this->genotypes = &genotypes->setImmutable();
    } else {
        this->genotypes = &GenoTypesContext::NO_GENOTYPES;
    }

    for(const std::shared_ptr<Allele> & allele : *alleles) {
        if(allele->getIsReference()) {
            this->REF = allele;
        } else if(alleles->size() == 2) {
            this->ALT = allele;
        }
    }

    this->fullyDecoded = fullyDecoded;
    if(!validationToPerform.empty()) {
        validate(validationToPerform);
    }
}

std::vector<std::shared_ptr<Allele>> VariantContext::makeAlleles(std::vector<std::shared_ptr<Allele>> &_alleles) {
    std::vector<std::shared_ptr<Allele>> alleleList;
    bool sawRef = false;
    for(const std::shared_ptr<Allele>& allele : _alleles) {
        int i = 0;
        for(int alleleListSize = _alleles.size(); i < alleleListSize; ++i) {
            if(i < alleleList.size() && alleleList.at(i) != nullptr && (*allele).equals(*alleleList.at(i), true)) {
                throw std::invalid_argument("Duplicate allele added to VariantContext");
            }
        }

        if(allele->getIsReference()) {
            if(sawRef) {
                throw std::invalid_argument("Alleles for a VariantContext must contain at most one reference allele");
            }
            alleleList.insert(alleleList.begin(), allele);
            sawRef = true;
        } else {
            alleleList.emplace_back(allele);
        }
    }

    if(alleleList.empty()) {
        throw std::invalid_argument("Cannot create a VariantContext with an empty allele list");
    } else if(alleleList.at(0)->getIsNonReference()) {
        throw std::invalid_argument("Alleles for a VariantContext must contain at least one reference allele");
    } else {
        return alleleList;
    }
}

bool VariantContext::validate(const std::set<Validation>& validationToPerform) {
    validateStop();
    for(Validation validation : validationToPerform) {
        switch (validation) {
            case ALLELES:
                validateAlleles();
                break;
            case GENOTYPES:
                validateGenotypes();
                break;
            default:
                throw std::invalid_argument("Unexpected validation mode ");
        }
    }
    return true;
}

bool VariantContext::hasAttribute(std::string &key) {
    return commonInfo.hasAttribute(key);
}

void VariantContext::validateStop() {
    if(hasAttribute((std::string &) "END")) {
        int end = getAttributeAsInt((std::string &) "END", -1);

        assert(end != -1);

        if(end != getEnd()) {
            throw std::invalid_argument("Badly formed variant context at location");
        }
    } else {
        long length = stop - start + 1;
        if(!hasSymbolicAlleles() && length != getReference()->getLength()){
            throw std::invalid_argument("BUG: GenomeLoc ");
        }
    }
}

int VariantContext::getAttributeAsInt(std::string &key, int defaultValue) {
    return commonInfo.getAttributeAsInt(key, defaultValue);
}

int VariantContext::getEnd() const {
    return (int)stop;
}

bool VariantContext::hasSymbolicAlleles() {
    return hasSymbolicAlleles(getAlleles());
}

std::vector<std::shared_ptr<Allele>> & VariantContext::getAlleles() {
    return alleles;
}

bool VariantContext::hasSymbolicAlleles(const std::vector<std::shared_ptr<Allele>> & alleles) {
    int i = 0;

    for(int size = (int)alleles.size(); i < size; i++) {
        if(alleles.at(i)->getIsSymbolic()) {
            return true;
        }
    }

    return false;
}

std::shared_ptr<Allele> VariantContext::getReference() {
    if(REF == nullptr) {
        throw std::invalid_argument("BUG: no reference allele found");
    } else {
        return REF;
    }
}

void VariantContext::validateAlleles() {
    bool alreadySeenRef = false;
    for(const std::shared_ptr<Allele> & allele : alleles) {
        if(allele->getIsReference()) {
            if(alreadySeenRef) {
                throw std::invalid_argument("BUG: Received two reference tagged alleles in VariantContext");
            }
            alreadySeenRef = true;
        }
        if(allele->getIsNoCall()) {
            throw std::invalid_argument("BUG: Cannot add a no call allele to a variant context");
        }
    }
    if(!alreadySeenRef) {
        throw std::invalid_argument("No reference allele found in VariantContext");
    }
}

void VariantContext::validateGenotypes() {
    if(genotypes == nullptr) {
        throw std::invalid_argument("Genotypes is null");
    } else {
        for(int i = 0; i < genotypes->getSize(); ++i) {
            Genotype* genotype = genotypes->get(i);
            if(genotype->isAvailable()) {
                std::vector<std::shared_ptr<Allele>> new_alleles = genotype->getAlleles();
                for(std::shared_ptr<Allele> allele : new_alleles) {
                    if(!hasAllele(allele) && allele->getIsCalled()) {
                        throw std::invalid_argument("Allele in genotype not in the variant context");
                    }
                }
            }
        }
    }
}

bool VariantContext::hasAllele(const std::shared_ptr<Allele> & allele) {
    return hasAllele(allele, false, true);
}

bool VariantContext::hasAllele(const std::shared_ptr<Allele> &allele, bool ignoreRefState) {
    return hasAllele(allele, ignoreRefState, true);
}

bool VariantContext::hasAllele(const std::shared_ptr<Allele> &allele, bool ignoreRefState, bool considerRefAllele) {
    if((!considerRefAllele || !((*allele) == (*REF))) && !((*allele) == (*ALT))) {
        std::vector<std::shared_ptr<Allele>> allelesToConsider = considerRefAllele ? getAlleles() : getAlternateAlleles();
        int i = 0;
        for(const std::shared_ptr<Allele> & allele1 : allelesToConsider) {
            if(allele1->equals(*allele, ignoreRefState)) {
                return true;
            }
        }
        return false;
    } else {
        return true;
    }
}

//TODO:验证写法是否正确
std::vector<std::shared_ptr<Allele>> VariantContext::getAlternateAlleles() {
    return {alleles.begin()++, alleles.end()};
}

int VariantContext::getStart() {
    return (int)start;
}

bool VariantContext::isBiallelic() {
    return getNAlleles() == 2;
}

int VariantContext::getNAlleles() {
    return (int)alleles.size();
}

bool VariantContext::isSNP() {
    return getType() == VariantContext_SNP;
}

VariantContextType VariantContext::getType() {
    if(type == VariantContext_NULL) {
        determineType();
    }
    return type;
}

void VariantContext::determineType() {
    if(type == VariantContext_NULL) {
        switch (getNAlleles()) {
            case 0:
                throw std::invalid_argument("Unexpected error: requested type of VariantContext with no alleles!");
            case 1:
                type = VariantContext_NO_VARIATION;
                break;
            default:
                determinePolymorphicType();
        }
    }
}

void VariantContext::determinePolymorphicType() {
    type = VariantContext_NULL;
    for(const std::shared_ptr<Allele>& allele : alleles) {
        if(allele != REF) {
            VariantContextType biallelicType = typeOfBiallelicVariant(REF, allele);
            if(type == VariantContext_NULL) {
                type = biallelicType;
            } else if (biallelicType != type) {
                type = VariantContext_MIXED;
                return;
            }
        }
    }
}

VariantContextType VariantContext::typeOfBiallelicVariant(const std::shared_ptr<Allele> & ref, const std::shared_ptr<Allele> &allele) {
    if(ref->getIsSymbolic()) {
        throw std::invalid_argument("Unexpected error: encountered a record with a symbolic reference allele");
    } else if (allele->getIsSymbolic()) {
        return VariantContext_SYMBOLIC;
    } else if (ref->getLength() == allele->getLength()) {
        return allele->getLength() == 1 ? VariantContext_SNP : VariantContext_MNP;
    } else {
        return VariantContext_INDEL;
    }
}

bool VariantContext::isSimpleDeletion() {
    return isSimpleIndel() && getAlternateAllele(0)->getLength() == 1;
}

std::shared_ptr<Allele> VariantContext::getAlternateAllele(int i) {
    return alleles.at(i+1);
}

bool VariantContext::isSimpleIndel() {
    return getType() == VariantContext_INDEL && isBiallelic() && getReference()->getLength() > 0 && getAlternateAllele(0)->getLength() > 0
    && getReference()->getBases().get()[0] == getAlternateAllele(0)->getBases().get()[0] && (getReference()->getLength() == 1 ||
            getAlternateAllele(0)->getLength() == 1);
}

bool VariantContext::isSimpleInsertion() {
    return isSimpleIndel() && getReference()->getLength() == 1;
}

std::map<std::string, void *> & VariantContext::getAttributes() {
    return commonInfo.getAttributes();
}

std::string &VariantContext::getContig() {
    return contig;
}

std::set<std::string> *VariantContext::getFiltersMaybeNull() {
    return commonInfo.getFiltersMaybeNull();
}

GenoTypesContext *VariantContext::getGenotypes() {
    return genotypes;
}

std::string &VariantContext::getID() {
    return ID;
}

double VariantContext::getLog10PError() {
    return commonInfo.getLog10PError();
}

std::string &VariantContext::getSource() {
    return commonInfo.getName();
}

bool VariantContext::isFullyDecoded() {
    return fullyDecoded;
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