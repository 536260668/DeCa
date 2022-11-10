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

    M2FiltersArgumentCollection();

private:
    static const double  DEFAULT_LOG_INDEL_PRIOR;
    static const double  DEFAULT_LOG_SNV_PRIOR;
    static const double  DEFAULT_LOG_INDEL_PRIOR_FOR_MITO;
    static const double  DEFAULT_LOG_SNV_PRIOR_FOR_MITO;
    static const double  DEFAULT_INITIAL_POSTERIOR_THRESHOLD;
    static const double  DEFAULT_F_SCORE_BETA;
    static const double  DEFAULT_MAX_FALSE_DISCOVERY_RATE;
};


#endif //MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H
