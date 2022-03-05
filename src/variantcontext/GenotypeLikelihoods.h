//
// Created by 梦想家xixi on 2021/11/27.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPELIKELIHOODS_H
#define MUTECT2CPP_MASTER_GENOTYPELIKELIHOODS_H

#include <string>
#include <map>
#include <utility>
#include "GenotypeLikelihoodsAllelePair.h"
#include "Genotype.h"
#include "VariantContext.h"
#include "GenotypeType.h"

class Genotype;

class GenotypeLikelihoods {
private:
    static const int NUM_LIKELIHOODS_CACHE_N_ALLELES = 5;
    static const int NUM_LIKELIHOODS_CACHE_PLOIDY = 10;
    static int numLikelihoodCache[5][10];
    double * log10Likelihoods = nullptr;
    int length;
    std::string likelihoodsAsString_PLs;
    static GenotypeLikelihoodsAllelePair** diploidPLIndexToAlleleIndex;
    static int allelePairLength;
    static double* parseDeprecatedGLString(std::string GLString, int &length);
    GenotypeLikelihoods(std::string asString) : likelihoodsAsString_PLs(std::move(asString)), length(0), log10Likelihoods(
            nullptr){}
    GenotypeLikelihoods(double * asVector, int length) : log10Likelihoods(asVector), length(length){}
    static double* PLsToGLs(int * pls, int length);
    static int * GLsToPLs(double* GLs, int length);
    static double* parsePLsIntoLikelihoods(std::string likelihoodsAsString_PLs, int & length);
    static double maxPL(double* GLs, int length);
    static std::string convertLikelihoodsToPLString(double* GLs, int length);
    double getLog10GQ(const std::vector<std::shared_ptr<Allele>> &genotypeAlleles, std::vector<Allele*> &contextAlleles);
    static GenotypeLikelihoodsAllelePair** calculateDiploidPLcache(int altAlleles, int & length);
    static int calcNumLikelihoods(int numAlleles, int ploidy);
    static void calculatePLIndexToAlleleIndices(int altAlleles, int ploidy, std::vector<std::vector<int>> &anyploidPLIndexToAlleleIndices, const std::vector<int>& genotype);

public:
    static const int MAX_PL = 2147483647;
    static const int MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED = 50;
    static GenotypeLikelihoods* fromPLField(std::string & PLs);
    static GenotypeLikelihoods* fromGLField(std::string & GLs);
    static GenotypeLikelihoods* fromLog10Likelihoods(double * log10Likelihoods, int length);
    static GenotypeLikelihoods* fromPLs(int* pls, int length);
    double* getAsVector();
    int* getAsPLs();
    std::string getAsString();
    std::map<GenotypeType, double> getAsMap(bool normalizeFromLog10);
    static double getGQLog10FromLikelihoods(int iOfChoosenGenotype, double * likelihoods, int length);
    double getLog10GQ(GenotypeType genotypeType);
    static int calculatePLindex(int allele1Index, int allele2Index);
    double getLog10GQ(Genotype* genotype, std::vector<Allele*> &contextAlleles);
    //TODO::double getLog10GQ(Genotype* genotype, VariantContext* context);
    static int numLikelihoods(int numAlleles, int ploidy);
    static GenotypeLikelihoodsAllelePair* getAllelePair(int PLindex);
    static void initializeAnyploidPLIndexToAlleleIndices(int altAlleles, int ploidy);
    static std::vector<int> getAlleles(int PLindex, int ploidy);
    static int* getPLIndicesOfAlleles(int allele1Index, int allele2Index);
    static void initial();

protected:
    static std::map<int, std::vector<std::vector<int>>> anyploidPloidyToPLIndexToAlleleIndices;
    static int* PLindexConversion;
    static std::vector<std::vector<int>> calculateAnyploidPLcache(int altAlleles, int ploidy);
};


#endif //MUTECT2CPP_MASTER_GENOTYPELIKELIHOODS_H
