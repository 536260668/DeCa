//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_VCFCONSTANTS_H
#define MUTECT2CPP_MASTER_VCFCONSTANTS_H

#include <string>

class VCFConstants {
public:
    // reserved INFO/FORMAT field keys
    inline const static std::string ANCESTRAL_ALLELE_KEY = "AA";    //---in c++17, you can use inline const static like Java
    inline const static std::string ALLELE_COUNT_KEY = "AC";
    inline const static std::string ALLELE_FREQUENCY_KEY = "AF";
    inline const static std::string ALLELE_NUMBER_KEY = "AN";
    inline const static std::string RMS_BASE_QUALITY_KEY = "BQ";
    inline const static std::string CIGAR_KEY = "CIGAR";
    inline const static std::string DBSNP_KEY = "DB";
    inline const static std::string DEPTH_KEY = "DP";
    inline const static std::string END_KEY = "END";

    inline const static std::string GENOTYPE_FILTER_KEY = "FT";
    inline const static std::string GENOTYPE_KEY = "GT";
    inline const static std::string GENOTYPE_POSTERIORS_KEY = "GP";
    inline const static std::string GENOTYPE_QUALITY_KEY = "GQ";
    inline const static std::string GENOTYPE_ALLELE_DEPTHS = "AD"; //AD isn't reserved, but is specifically handled by VariantContext
    inline const static std::string GENOTYPE_PL_KEY = "PL";   // phred-scaled genotype likelihoods
    inline const static std::string EXPECTED_ALLELE_COUNT_KEY = "EC";


    // missing/default values
    inline const static std::string UNFILTERED = ".";
    inline const static std::string PASSES_FILTERS_v3 = "0";
    inline const static std::string PASSES_FILTERS_v4 = "PASS";
    inline const static std::string EMPTY_ID_FIELD = ".";
    inline const static std::string EMPTY_INFO_FIELD = ".";
    inline const static std::string EMPTY_ALTERNATE_ALLELE_FIELD = ".";
    inline const static std::string MISSING_VALUE_v4 = ".";
    inline const static std::string MISSING_QUALITY_v3 = "-1";
    inline const static double MISSING_QUALITY_v3_DOUBLE = stod(MISSING_QUALITY_v3);
};


#endif //MUTECT2CPP_MASTER_VCFCONSTANTS_H
