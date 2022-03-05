//
// Created by 梦想家xixi on 2021/11/27.
//

#include "GenotypeLikelihoods.h"
#include <regex>
#include <utility>
#include <vector>
#include <cmath>
#include <sstream>
#include "StringUtils.h"
#include "utils/GeneralUtils.h"

int GenotypeLikelihoods::numLikelihoodCache[5][10] = {0};

int GenotypeLikelihoods::allelePairLength = 0;

GenotypeLikelihoodsAllelePair** GenotypeLikelihoods::diploidPLIndexToAlleleIndex = nullptr;

std::map<int, std::vector<std::vector<int>>> GenotypeLikelihoods::anyploidPloidyToPLIndexToAlleleIndices;

int* GenotypeLikelihoods::PLindexConversion = nullptr;

GenotypeLikelihoods *GenotypeLikelihoods::fromGLField(std::string & GLs) {
    int length;
    double * vec = parseDeprecatedGLString(GLs, length);
    return new GenotypeLikelihoods(vec, length);
}

double *GenotypeLikelihoods::parseDeprecatedGLString(std::string GLString, int& length) {
    if(GLString == ".")
        return nullptr;
    else {
        std::regex ws_re(",");
        std::vector<std::string> strings(std::sregex_token_iterator(GLString.begin(), GLString.end(), ws_re, -1), std::sregex_token_iterator());
        length = strings.size();
        double * likelihoodsAsVector = new double[length];
        int i = 0;
        for(std::string str : strings) {
            likelihoodsAsVector[i] = StringUtils::parseDouble(str);
        }
        return likelihoodsAsVector;
    }
}

GenotypeLikelihoods *GenotypeLikelihoods::fromLog10Likelihoods(double *log10Likelihoods, int length) {
    return new GenotypeLikelihoods(log10Likelihoods, length);
}

GenotypeLikelihoods *GenotypeLikelihoods::fromPLField(std::string &PLs) {
    return new GenotypeLikelihoods(PLs);
}

double *GenotypeLikelihoods::PLsToGLs(int *pls, int length) {
    double * likelihoodsAsVector = new double[length];
    for(int i = 0; i < length; i++) {
        likelihoodsAsVector[i] = (double)pls[i] / -10.0;
    }
    return likelihoodsAsVector;
}

GenotypeLikelihoods *GenotypeLikelihoods::fromPLs(int *pls, int length) {
    return new GenotypeLikelihoods(PLsToGLs(pls, length), length);
}

double *GenotypeLikelihoods::getAsVector() {
    if(log10Likelihoods == nullptr) {
        log10Likelihoods = parsePLsIntoLikelihoods(likelihoodsAsString_PLs, length);
    }
    return log10Likelihoods;
}

double *GenotypeLikelihoods::parsePLsIntoLikelihoods(std::string likelihoodsAsString_PLs, int & length) {
    if(likelihoodsAsString_PLs != ".") {
        std::regex ws_re(",");
        std::vector<std::string> strings(std::sregex_token_iterator(likelihoodsAsString_PLs.begin(), likelihoodsAsString_PLs.end(), ws_re, -1), std::sregex_token_iterator());
        length = strings.size();
        double* likelihoodsAsVector = new double[length];
        for(int i = 0; i < length; i++) {
            likelihoodsAsVector[i] = StringUtils::parseInt(strings[i]) / -10.0;
        }
        return likelihoodsAsVector;
    } else
        return nullptr;
}

int *GenotypeLikelihoods::getAsPLs() {
    double * GLs = getAsVector();
    return GLs == nullptr ? nullptr : GLsToPLs(GLs, length);
}

int *GenotypeLikelihoods::GLsToPLs(double *GLs, int length) {
    int * pls = new int[length];
    double adjust = maxPL(GLs, length);
    for(int i = 0; i < length; ++i) {
        pls[i] = (int)std::round(std::min(-10.0 * (GLs[i] - adjust), 2.147483647E9));
    }
    return pls;
}

double GenotypeLikelihoods::maxPL(double *GLs, int length) {
    double adjust = -1.0 / 0.0;
    double * var3 = GLs;
    int var4 = length;

    for(int var5 = 0 ; var5 < var4; ++var5) {
        double l = var3[var5];
        adjust = std::max(adjust, l);
    }
    return adjust;
}

std::string GenotypeLikelihoods::getAsString() {
    if(likelihoodsAsString_PLs.empty()) {
        if(log10Likelihoods == nullptr) {
            throw std::invalid_argument("BUG: Attempted to get likelihoods as strings and neither the vector nor the string is set!");
        }
        likelihoodsAsString_PLs = convertLikelihoodsToPLString(log10Likelihoods, length);
    }

    return likelihoodsAsString_PLs;
}

std::string GenotypeLikelihoods::convertLikelihoodsToPLString(double *GLs, int length) {
    if(GLs == nullptr) {
        return ".";
    } else {
        std::stringstream ss;
        ss.precision(10);
        bool first = true;
        int* var3 = GLsToPLs(GLs, length);
        int var4 = length;

        for(int var5 = 0; var5 < var4; var5++) {
            int pl = var3[var5];
            if(!first) {
                ss << ',';
            } else {
                first  = false;
            }
            ss << pl;
        }
        std::string ret;
        ss >> ret;
        return ret;
    }
}

std::map<GenotypeType, double> GenotypeLikelihoods::getAsMap(bool normalizeFromLog10) {
    double * asVector = getAsVector();
    double* likelihoods = normalizeFromLog10 ? GeneralUtils::normalizeFromLog10(asVector, length) : asVector;
    if(likelihoods == nullptr) {
        return {};
    } else {
        std::map<GenotypeType, double> likelihoodsMap;
        likelihoodsMap.insert(std::pair<GenotypeType, double>(HOM_REF, likelihoods[HOM_REF-1]));
        likelihoodsMap.insert(std::pair<GenotypeType, double>(HET, likelihoods[HET-1]));
        likelihoodsMap.insert(std::pair<GenotypeType, double>(HOM_VAR, likelihoods[HOM_VAR-1]));
        return likelihoodsMap;
    }
}

double GenotypeLikelihoods::getGQLog10FromLikelihoods(int iOfChoosenGenotype, double *likelihoods, int length) {
    if(likelihoods == nullptr) {
        return -1.0 / 0.0;
    } else {
        double qual = -1.0 / 0.0;

        for(int i = 0; i < length; ++i) {
            if(i != iOfChoosenGenotype && likelihoods[i] >= qual) {
                qual = likelihoods[i];
            }
        }
        qual = likelihoods[iOfChoosenGenotype] - qual;
        if(qual < 0.0) {
            double * normalized = GeneralUtils::normalizeFromLog10(likelihoods, length);
            double chosenGenotype = normalized[iOfChoosenGenotype];
            return std::log10(1.0 - chosenGenotype);
        } else {
            return -1.0 * qual;
        }
    }
}

double GenotypeLikelihoods::getLog10GQ(GenotypeType genotypeType) {
    double * asVector = getAsVector();
    return getGQLog10FromLikelihoods(genotypeType - 1, asVector, length);
}

double GenotypeLikelihoods::getLog10GQ(const std::vector<std::shared_ptr<Allele>> &genotypeAlleles, std::vector<Allele *> &contextAlleles) {
    int allele1Index = 0;
    int allele2Index = 0;
    int i = 0;
    for(Allele* allele : contextAlleles) {
        if((*genotypeAlleles.at(0)) == (*allele)) {
            allele1Index = i;
            break;
        }
        i++;
    }
    for(Allele* allele : contextAlleles) {
        if((*genotypeAlleles.at(1)) == (*allele)) {
            allele2Index = i;
            break;
        }
        i++;
    }
    int plIndex = calculatePLindex(allele1Index, allele2Index);
    double * asVector = getAsVector();
    return getGQLog10FromLikelihoods(plIndex, getAsVector(), length);
}

int GenotypeLikelihoods::calculatePLindex(int allele1Index, int allele2Index) {
    return allele2Index * (allele2Index + 1) / 2 + allele1Index;
}

double GenotypeLikelihoods::getLog10GQ(Genotype *genotype, std::vector<Allele *> &contextAlleles) {
    std::vector<std::shared_ptr<Allele>> genotypes = genotype->getAlleles();
    return getLog10GQ(genotypes, contextAlleles);
}

GenotypeLikelihoodsAllelePair** GenotypeLikelihoods::calculateDiploidPLcache(int altAlleles, int & length) {
    int numLikelihood = numLikelihoods(1 + altAlleles, 2);
    int i;
    GenotypeLikelihoodsAllelePair** cache = new GenotypeLikelihoodsAllelePair*[numLikelihood];
    for(i = 0; i <= altAlleles; ++i) {
        for(int allele2 = i; allele2 <= altAlleles; ++allele2) {
            cache[calculatePLindex(i, allele2)] = new GenotypeLikelihoodsAllelePair(i, allele2);
        }
    }

    for(i = 0; i < numLikelihood; ++i) {
        if(cache[i] == nullptr)
            throw std::invalid_argument("BUG: cache entry is unexpected null");
    }
    length = numLikelihood;
    return cache;

}

int GenotypeLikelihoods::numLikelihoods(int numAlleles, int ploidy) {
    return numAlleles < 5 && ploidy < 10 ? numLikelihoodCache[numAlleles][ploidy] : calcNumLikelihoods(numAlleles, ploidy);
}

int GenotypeLikelihoods::calcNumLikelihoods(int numAlleles, int ploidy) {
    if (numAlleles == 1) {
        return 1;
    } else if (ploidy == 1) {
        return numAlleles;
    } else {
        int acc = 0;

        for(int k = 0; k <= ploidy; ++k) {
            acc += calcNumLikelihoods(numAlleles - 1, ploidy - k);
        }

        return acc;
    }
}

void GenotypeLikelihoods::calculatePLIndexToAlleleIndices(int altAlleles, int ploidy,
                                                          std::vector<std::vector<int>> &anyploidPLIndexToAlleleIndices,
                                                          const std::vector<int>& genotype) {
    for(int a = 0; a <= altAlleles; ++a) {
        std::vector<int> gt;
        gt.emplace_back(a);
        for(int tmp : genotype) {gt.emplace_back(tmp);}
        if (ploidy == 1) {
            anyploidPLIndexToAlleleIndices.emplace_back(gt);
        } else if (ploidy > 1) {
            calculatePLIndexToAlleleIndices(a, ploidy - 1, anyploidPLIndexToAlleleIndices, gt);
        }
    }
}

std::vector<std::vector<int>> GenotypeLikelihoods::calculateAnyploidPLcache(int altAlleles, int ploidy) {
    std::vector<std::vector<int>> anyploidPLIndexToAlleleIndices;
    calculatePLIndexToAlleleIndices(altAlleles, ploidy, anyploidPLIndexToAlleleIndices, std::vector<int>());
    return anyploidPLIndexToAlleleIndices;
}

GenotypeLikelihoodsAllelePair *GenotypeLikelihoods::getAllelePair(int PLindex) {
    if(PLindex >= 0 && PLindex < allelePairLength) {
        return diploidPLIndexToAlleleIndex[PLindex];
    } else {
        throw std::invalid_argument("The PL index cannot be negative");
    }
}

//TODO:需要同步
void GenotypeLikelihoods::initializeAnyploidPLIndexToAlleleIndices(int altAlleles, int ploidy) {
    if(altAlleles <= 0) {
        throw std::invalid_argument("Must have at least one alternate allele");
    } else if (ploidy <= 0) {
        throw std::invalid_argument("Ploidy must be at least 1");
    } else {
        anyploidPloidyToPLIndexToAlleleIndices.insert(std::pair<int, std::vector<std::vector<int>>>(ploidy, calculateAnyploidPLcache(altAlleles, ploidy)));
    }
}

//TODO:需要同步
std::vector<int> GenotypeLikelihoods::getAlleles(int PLindex, int ploidy) {
    if(ploidy == 2) {
        GenotypeLikelihoodsAllelePair* tmpPair = getAllelePair(PLindex);
        return std::vector<int>{tmpPair->alleleIndex1, tmpPair->alleleIndex2};
    } else if (anyploidPloidyToPLIndexToAlleleIndices.find(ploidy) == anyploidPloidyToPLIndexToAlleleIndices.end()) {
        throw std::invalid_argument("Must initialize the cache of allele anyploid indices for ploidy ");
    } else if (PLindex >= 0 && PLindex < anyploidPloidyToPLIndexToAlleleIndices.at(ploidy).size()) {
        return anyploidPloidyToPLIndexToAlleleIndices.at(ploidy).at(PLindex);
    } else {
        throw std::invalid_argument("cannot have a negative value.");
    }
}

int *GenotypeLikelihoods::getPLIndicesOfAlleles(int allele1Index, int allele2Index) {
    int* indexes = new int[]{calculatePLindex(allele1Index, allele1Index), calculatePLindex(allele1Index, allele2Index), calculatePLindex(allele2Index, allele2Index)};
    return indexes;
}

void GenotypeLikelihoods::initial() {
    for(int numAlleles = 1; numAlleles < 5; ++numAlleles) {
        for(int ploidy = 1; ploidy < 10; ++ploidy) {
            numLikelihoodCache[numAlleles][ploidy] = calcNumLikelihoods(numAlleles, ploidy);
        }
    }

    diploidPLIndexToAlleleIndex = calculateDiploidPLcache(50, allelePairLength);
    PLindexConversion = new int[]{0, 1, 3, 6, 2, 4, 7, 5, 8, 9};
}
