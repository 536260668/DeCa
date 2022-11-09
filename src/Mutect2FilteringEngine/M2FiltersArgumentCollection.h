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

private:
    static const double  DEFAULT_LOG_INDEL_PRIOR;
    static const double  DEFAULT_LOG_SNV_PRIOR;
    static const double  DEFAULT_LOG_INDEL_PRIOR_FOR_MITO;
    static const double  DEFAULT_LOG_SNV_PRIOR_FOR_MITO;
    M2FiltersArgumentCollection();
};


#endif //MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H
