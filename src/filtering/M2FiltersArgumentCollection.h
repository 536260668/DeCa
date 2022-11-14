//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H


class M2FiltersArgumentCollection {
public:
    double getLogIndelPrior();
    double getLogSnvPrior();
    bool mitochondria;
    double logIndelPrior;
    double logSNVPrior;
    double initialPosteriorThreshold;
    double fScoreBeta;
    double maxFalsePositiveRate;
    double initialLogPriorOfVariantVersusArtifact;
    int minMedianBaseQuality;
    int minMedianMappingQuality = DEFAULT_MIN_MEDIAN_MAPPING_QUALITY;
    int longIndelLength = DEFAULT_LONG_INDEL_SIZE;

    M2FiltersArgumentCollection();

private:
    static const double  DEFAULT_LOG_INDEL_PRIOR;
    static const double  DEFAULT_LOG_SNV_PRIOR;
    static const double  DEFAULT_LOG_INDEL_PRIOR_FOR_MITO;
    static const double  DEFAULT_LOG_SNV_PRIOR_FOR_MITO;
    static const double  DEFAULT_INITIAL_POSTERIOR_THRESHOLD;
    static const double  DEFAULT_F_SCORE_BETA;
    static const double  DEFAULT_MAX_FALSE_DISCOVERY_RATE;
    static const double  DEFAULT_INITIAL_LOG_PRIOR_OF_VARIANT_VERSUS_ARTIFACT;
    static const int DEFAULT_MIN_MEDIAN_BASE_QUALITY = 20;
    static const int DEFAULT_MIN_MEDIAN_MAPPING_QUALITY = 30;
    static const int DEFAULT_LONG_INDEL_SIZE = 5;
};


#endif //MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H
