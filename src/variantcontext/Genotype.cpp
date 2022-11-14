//
// Created by lhh on 11/12/21.
//

#include "Genotype.h"
#include "Mutect2Utils.h"
#include <utility>
#include <stdexcept>
#include <algorithm>
#include "StringUtils.h"

const std::string Genotype::PRIMARY_KEYS[6] = {"FT", "GT", "GQ", "DP", "AD", "PL"};

const std::string Genotype::PHASED_ALLELE_SEPARATOR = "|";

const std::string Genotype::UNPHASED_ALLELE_SEPARATOR = "/";

Genotype::Genotype(std::string sampleName, const std::string& filters) {
    this->sampleName = std::move(sampleName);
    this->filters = !filters.empty() ? filters : "";
}

int Genotype::countAllele(const std::shared_ptr<Allele> &allele) {
    int c = 0;
    for(const std::shared_ptr<Allele> & var : getAlleles()) {
        if((*var) == (*allele))
            ++c;
    }
    return c;
}

GenotypeType Genotype::getType() {
    if(type == GenotypeType_NULL) {
        type = determineType();
    }
    return type;
}

GenotypeType Genotype::determineType() {
    std::vector<std::shared_ptr<Allele>> alleles = getAlleles();
    if(alleles.empty()) {
        return UNAVAILABLE;
    } else {
        bool sawNoCall = false;
        bool sawMultipleAlleles = false;
        std::shared_ptr<Allele> firstCallAllele = nullptr;

        for(const std::shared_ptr<Allele> & allele : alleles) {
            if(allele->getIsNoCall()) {
                sawNoCall = true;
            } else if (firstCallAllele == nullptr) {
                firstCallAllele = allele;
            } else if (allele != firstCallAllele) {
                sawMultipleAlleles = true;
            }
        }

        if(sawNoCall) {
            if(firstCallAllele == nullptr) {
                return NO_CALL;
            } else {
                return MIXED;
            }
        } else if (firstCallAllele == nullptr) {
            throw std::invalid_argument("BUG: there are no alleles present in this genotype but the alleles list is not null");
        } else {
            return sawMultipleAlleles ? HET : (firstCallAllele->getIsReference() ? HOM_REF : HOM_VAR);
        }
    }
}

bool Genotype::isHetNonRef() {
    return getType() == HET && getAllele(0)->getIsNonReference() && getAllele(1)->getIsNonReference();
}

GenotypeLikelihoods *Genotype::getLikelihoods() {
    int length;
    return hasLikelihoods() ? GenotypeLikelihoods::fromPLs(getPL(length), length) : nullptr;
}

bool Genotype::isNonInformative() {
    int length;
    if(getPL(length) == nullptr) {
        return true;
    } else {
        int* var1 = getPL(length);
        int var2 = length;

        for(int var3 = 0; var3 < var2; ++var3) {
            int PL = var1[var3];
            if (PL != 0) {
                return false;
            }
        }

        return true;
    }
}

bool Genotype::hasPL() {
    int length;
    return getPL(length) != nullptr;
}

bool Genotype::hasAD() {
    int length;
    return getAD(length) != nullptr;
}

bool Genotype::hasLikelihoods() {
    int length;
    return getPL(length) != nullptr;
}

std::string Genotype::getGenotypeString(bool ignoreRefState) {
    if(getPloidy() == 0) {
        return "NA";
    } else {
        std::string separator = isPhased() ? "|" : "/";
        if(ignoreRefState) {
            std::vector<std::string> lists = getAlleleStrings();
            return StringUtils::join(separator, lists);
        } else {
            std::vector<std::shared_ptr<Allele>> alleles = getAlleles();
            if(!isPhased()) {
                std::sort(alleles.begin(), alleles.end(), [&](const std::shared_ptr<Allele> & a1, const std::shared_ptr<Allele> & a2)->bool {return (*a1) < (*a2);});
            }
            std::vector<std::string> lists(alleles.size());
            for(const std::shared_ptr<Allele>& allele : alleles) {
                lists.emplace_back(allele->getBaseString());
            }
            return StringUtils::join(separator, lists);
        }
    }
}

std::vector<std::string> Genotype::getAlleleStrings() {
    std::vector<std::string> al(getPloidy());
    for(const std::shared_ptr<Allele>& allele : getAlleles()) {
        al.emplace_back(allele->getBaseString());
    }
    return al;
}

bool Genotype::operator<(const Genotype &other) const {
    return getSampleName() < other.getSampleName();
}

bool Genotype::sameGenotype(Genotype *other, bool ignorePhase) {
    if(getPloidy() != other->getPloidy()) {
        return false;
    } else {
        std::vector<std::shared_ptr<Allele>> thisAlleles = getAlleles();
        std::vector<std::shared_ptr<Allele>> otherAlleles = other->getAlleles();
        if(ignorePhase) {
            std::set<std::shared_ptr<Allele>> treeThisAlleles(thisAlleles.begin(), thisAlleles.end());
            std::set<std::shared_ptr<Allele>> treeOtherAlleles(otherAlleles.begin(), otherAlleles.end());
            return Mutect2Utils::isSetEquals(treeOtherAlleles, treeThisAlleles);
        }
        return Mutect2Utils::isVectorEquals(thisAlleles, otherAlleles);
    }
}

bool Genotype::sameGenotype(Genotype *other) {
    return sameGenotype(other, true);
}

bool Genotype::hasExtendedAttribute(const std::string &key) {
    std::map<std::string, AttributeValue> tmp = getExtendedAttributes();
    return tmp.find(key) != tmp.end();
}

AttributeValue Genotype::getExtendedAttribute(const std::string& key, void *defaultValue) {
    return hasExtendedAttribute(key) ? getExtendedAttributes().at(key) : *(AttributeValue *)defaultValue;
}

AttributeValue Genotype::getExtendedAttribute(const std::string & key) {
    return getExtendedAttribute(key, nullptr);
}

//需要delete
AttributeValue Genotype::getAnyAttribute(const std::string& key) {
    if(key == "GT") {
        std::vector<std::shared_ptr<Allele>> ret = getAlleles();
        return (new std::vector<std::shared_ptr<Allele>>(ret.begin(), ret.end()));
    } else if (key == "GQ") {
        int* ret = new int(getGQ());
        return ret;
    }  else {
        int* var3;
        int var4;
        int var5;
        int i;
        if(key == "AD") {
            if(!hasAD()) {
                return AttributeValue::empty_value();
            } else {
                int length;
                var3 = getAD(length);
                std::vector<int>* ret = new std::vector<int>(length);
                var4 = length;
                for(var5 = 0; var5 < var4; ++var5) {
                    i = var3[var5];
                    ret->emplace_back(i);
                }
                return ret;
            }
        } else if (key != "PL") {
            if(key == "DP") {
                return new int(getDP());
            } else {
                return key == "FT" ? new std::string(getFilters()) : getExtendedAttribute(key);
            }
        }else if (!hasPL()) {
            return AttributeValue::empty_value();;
        } else {
            int length;
            var3 = getPL(length);
            std::vector<int>* ret = new std::vector<int>(length);
            var4 = length;
            for(var5 = 0; var5 < var4; ++var5) {
                i = var3[var5];
                ret->emplace_back(i);
            }
            return ret;
        }
    }
}

bool Genotype::hasAnyAttribute(std::string key) {
    if(key == "GT") {
        return isAvailable();
    } else if (key == "GQ") {
        return hasGQ();
    } else if (key == "AD") {
        return hasAD();
    } else if (key == "PL") {
        return hasPL();
    } else if (key == "DP") {
        return hasDP();
    } else {
        return key == "FT" || hasExtendedAttribute(key);
    }
}

Genotype::~Genotype() = default;


